#=
    NAFF
    Julia v1.0+
=#
module NAFF

using FFTW
export naff

"""
naff(data, [turns=300, nterms=1, skipTurns=0, window=1])


# Description
Function that takes an Array (Complex or Float) and runs the NAFF [J. Laskar]. Implementation of the PyNAFF code [https://github.com/nkarast/PyNAFF](https://github.com/nkarast/PyNAFF)
algorithm.

# Usage & Arguments
- `data::Array{Float64|ComplexF64}` : The data on which to perform NAFF
- `turns::Int64` : Number of turns (i.e. indices) to use from data
- `nterms::Int64`: Number of frequencies to return. The code can return up to the maximum number of frequencies found in the signal, which can be less than nterms provided.
- `skipTurns::Int64` : Number of turns to skip. Data array starts from skipTurns+1
- `window::Int64` : Power of Hann window to use.

Returns a matrix. Each colum includes 5 rows of coded elements as
- order of harmonic
- frequency
- Amplitude
- Re{Amplitude}
- Im{Amplitude}

# Example:
```@example
julia> data = sin.(collect(0:999)*2*pi*0.123) .+ sin.(collect(0:999)*2*pi*0.234);

julia> NAFF.naff(data, 300, 4, 0, 1)
5×4 Array{Float64,2}:
  1.0           2.0           3.0         4.0
  0.123001     -0.123001      0.234      -0.234
  0.5           0.5           0.5         0.5
 -0.000431821  -0.000434902  -0.0003409  -1.45175e-7
 -0.5           0.5          -0.5         0.5
```

# Contact
N. Karastathis, nkarast .at. cern .dot. ch

# References:
J. Laskar, "The chaotic motion of the solar system: A numerical estimate of the size of the chaotic zones", Icarus Vol. 88, 2, Dec.1990, p.266-291
"""
function naff(data, turns=300, nterms=1, skipTurns=0, window=1)

    if turns >= length(data)
        error("#naff: Input data must be at least of length turns+1.")
    end

    if turns < 6
        error("#naff: Minimum number of turns is 6.")
    end

    if mod(turns,6) != 0
        a,b = divrem(turns,6)
        turns = Int(6*a)
    end

    NFR   = 1000
	vars  = Dict(
	"NFS"   => 0,
    "TFS"   => zeros(Float64, NFR),
    "ZAMP"  => zeros(ComplexF64, NFR),
    "ZALP"  => zeros(ComplexF64, (NFR,NFR)),
    "ZTABS" => ComplexF64[],
    "TWIN"  => Float64[],
	)

    function getIntegral(FR, turns)
        if mod(turns,6)!=0
            error("Turns need to be *6")
        end

        K = Int64(round((turns-1)/6))
        i_line = collect(1:turns)
        ZTF_tmp = vars["ZTABS"][2:end].*vars["TWIN"][2:end].*exp.(-2.0*pi*1.0im*FR .* i_line)
        ZTF = ComplexF64[vars["ZTABS"][1]*vars["TWIN"][1]]
        ZTF = append!(ZTF, ZTF_tmp)
        let ZOM = 41.0 * ZTF[1]+216.0* ZTF[2]+27.0* ZTF[3]+272.0* ZTF[4]+27.0* ZTF[5]+216.0*ZTF[6]+41.0*ZTF[end]
        for I in 1:(K-1)
            ZOM=ZOM+82.0*ZTF[6*I+1]+216.0*ZTF[6*I+2]+27.0*ZTF[6*I+3]+272.0*ZTF[6*I+4]+27.0*ZTF[6*I+5]+216.0*ZTF[6*I+6]
        end
        ZOM=ZOM*(1.0/turns)*(6.0/840.0)
        A = real(ZOM)
        B = imag(ZOM)
        RMD = abs(ZOM)
        return RMD, A, B
		end
    end #getIntegral

    function frefin(turns, FR, STAREP, EPS)
        EPSI = 1.0e-15
        X2  = FR
        PAS = STAREP
        Y2, A2, B2  = getIntegral(X2, turns)
        X1  = X2 - PAS
        X3  = X2 + PAS
        Y1, A1, B1  = getIntegral(X1, turns)
        Y3, A3, B3  = getIntegral(X3, turns)
        while true
        	if PAS >=EPS
        		if abs(Y3-Y1) < EPSI
        			break
                end
        		if (Y1<Y2) && (Y3<Y2)
        			R2  = (Y1-Y2)/(X1-X2)
        			R3  = (Y1-Y3)/(X1-X3)
        			A   = (R2 - R3)/(X2-X3)
        			B   = R2 - A*(X1+X2)
        			XX2 = -B/(2.0*A)
        			PAS = abs(XX2-X2)
        			if XX2 > X2
        				X1 = X2
        				Y1, A1, B1 = Y2, A2, B2
        				X2 = XX2
        				Y2, A2, B2 = getIntegral(X2, turns)
        				X3 = X2 + PAS
        				Y3, A3, B3 = getIntegral(X3, turns)
        			else
        				X3 = X2
        				Y3, A3, B3 = Y2, A2, B2
        				X2 = XX2
        				Y2, A2, B2 = getIntegral(X2, turns)
        				X1 = X2 - PAS
        				Y1, A1, B1 = getIntegral(X1, turns)
                    end
        		else
        			if Y1>Y3
        				X2 = X1
        				Y2, A2, B2 = Y1, A1, B1
        			else
        				X2 = X3
        				Y2, A2, B2 = Y3, A3, B3
                    end

        			X1 = X2 - PAS
        			X3 = X2 + PAS
        			Y1, A1, B1 = getIntegral(X1, turns)
        			Y3, A2, B2 = getIntegral(X3, turns)
        			if (Y3-Y1)-(Y3-Y2)==0.0
        				PAS=PAS+EPS
                    end
                end

        	else
        		break
            end
        end
        return X2, Y2, A2, B2
    end # frefin

    function fretes(FR, FREFON)
		#=
		If more than one term found, check how different they are
		=#
		TOL   = 1.0e-4 # this is defined in mftnaf in lashkar
		IFLAG = 1
		NUMFR = 0
		ECART = abs(FREFON)
		for i in 1:length(vars["TFS"])
			TEST = abs(vars["TFS"][i] - FR)
			if TEST < ECART
				if Float64(TEST)/Float64(ECART) < TOL
					IFLAG = -1
					NUMFR = i
					break
				else
					IFLAG = 0
					continue
                end
            end
        end #for
		return IFLAG, NUMFR
    end #fretes


    function modfre(turns, FR, NUMFR, A, B)
		#=
		If I found something very close to one of the FR before, I assume that this comes from data
		I had not removed successfully => Remove them without orthonormalization
		=#
		ZI  = 0. + 1.0im
		ZOM = 1.0im*FR
		ZA  = 1.0*A + 1.0im*B
		if length(vars["ZAMP"])<= NUMFR
			vars["ZAMP"][NUMFR] = ComplexF64(0.)
        end
		vars["ZAMP"][NUMFR] = vars["ZAMP"][NUMFR] .+ ZA
		i_line = collect(1, turns)
		ZT_tmp = ZA*exp.(2.0*(i_line)*pi*ZOM)
		ZT     = ComplexF64[ZA]
		ZT     = append!(ZT, ZT_tmp)
		ZTABS_tmp = vars["ZTABS"] .- ZT
		vars["ZTABS"] = ZTABS_tmp
    end #modfre

    function proscaa(turns, FS, FS_OLD)  # ok
		ZI = 0.0+1.0im
		OM = FS-FS_OLD
		ANGI = 2.0*pi*OM
		i_line = collect(1:turns)
		ZT_tmp = exp.(-2.0*(i_line)*1.0im*pi*OM)
		ZT_zero = ComplexF64[1.0]
		ZT = append!(ZT_zero, ZT_tmp)
		ZTF = vars["TWIN"] .* ZT
		let ZOM = 41.0*ZTF[1]+216.0*ZTF[2]+27.0*ZTF[3]+272.0*ZTF[4]+27.0*ZTF[5]+216.0*ZTF[6]+41.0*ZTF[turns]
		for I in 1:Int64(round((turns-1)/6))-1
			ZOM=ZOM+82.0*ZTF[6*I+1]+216.0*ZTF[6*I+2]+27.0*ZTF[6*I+3]+272.0*ZTF[6*I+4]+27.0*ZTF[6*I+5]+216.0*ZTF[6*I+6]
        end

		ZOM=ZOM*(1.0/turns)*(6.0/840.0)
		return ZOM
		end # let
    end # proscaa

    function gramsc(turns, FR, A, B)
		#=
		Remove the contribution of the frequency found from the Data and orthonormalize
		=#
		ZTEE = zeros(ComplexF64, vars["NFS"]+1)
		for i in 1:vars["NFS"]
			ZTEE[i] = proscaa(turns, FR, vars["TFS"][i])
        end
		NF = vars["NFS"]+1
		ZTEE[NF] = 1.0+0.0im
		vars["TFS"][NF] = FR
		for k in 1:vars["NFS"]
			for i in k:vars["NFS"]
				for j in 1:i
					vars["ZALP"][NF, k] = vars["ZALP"][NF, k] - conj(vars["ZALP"][i,j])*vars["ZALP"][i,k]*ZTEE[j]
                end
            end
        end

		vars["ZALP"][NF, NF] = 1.0+0.0im
		DIV  = 1.0
		ZDIV = 0.0+0.0im
		for i in 1:NF
			ZDIV = ZDIV + conj(vars["ZALP"][NF, i])*ZTEE[i]
        end
		DIV = sqrt(abs(ZDIV))
		vars["ZALP"][NF,:] = vars["ZALP"][NF,:]/DIV
		ZMUL = ComplexF64(A,B)/DIV
		ZI = 0.0+1.0im

        i_line = collect(1:turns)
		for i in 1:NF
			ZOM = 1.0im*vars["TFS"][i]
			ZA  = vars["ZALP"][NF,i]*ZMUL
			vars["ZAMP"][i] = vars["ZAMP"][i] .+ ZA
			ZT_zero = ComplexF64[ZA]
			ZT_tmp = ZA*exp.(2.0*(i_line)*1.0im*pi*vars["TFS"][i])
			ZT = append!(ZT_zero, ZT_tmp)
			vars["ZTABS"] = vars["ZTABS"] .- ZT
        end
    end #gramsc
    # -------------------------

    FREFON = 1.0/turns
    NEPS = 100000000
    EPS = FREFON/NEPS

    T = collect(0:turns)*2.0*pi .- pi*turns
    vars["TWIN"]  = 1.0.+cos.(T/turns)
	vars["TWIN"]  = ((2.0^window*factorial(window)^2)/Float64(factorial(2*window)))*(1.0.+cos.(T/turns)).^window
	vars["ZTABS"] = data[skipTurns+1:skipTurns+turns+1] # ztabs[end] = data[turns]

    TOL = 1.0e-4
	STAREP = FREFON/3.0
	for term in 1:nterms
		data_for_fft = vars["ZTABS"] .* vars["TWIN"]
        data_for_fft = data_for_fft[1:length(data_for_fft)-1]

		y = fft(data_for_fft)

		RTAB = sqrt.(real(y).^2 + imag(y).^2)/turns  # normalized
		INDX = argmax(RTAB)
		VMAX = maximum(RTAB)

		if INDX == 1
			println("## naff: Remove the DC component from the data (i.e. the mean).")
        end
		if INDX <= turns/2.0
			IFR = INDX - 1
		else
			IFR = INDX-1-turns
        end

		FR = (IFR)*FREFON
		FR, RMD, A, B = frefin(turns, FR, STAREP, EPS)
		IFLAG, NUMFR = fretes(FR, FREFON)
		if IFLAG == 1
			gramsc(turns, FR, A, B)
			vars["NFS"] = vars["NFS"] + 1
		elseif IFLAG == 0
			# continue
			break  # if I put continue it will find again and again the same freq/ with break it stops repeating
		elseif IFLAG == -1
			modfre(turns, FR, NUMFR, A, B)
        end
    end


	let result = Float64[]

	for i in 1:vars["NFS"]
		if i==1
			AMP = abs(vars["ZAMP"][1])
			result = vcat(result, [Int(1), vars["TFS"][1], AMP, real(vars["ZAMP"][1]), imag(vars["ZAMP"][1])])
		else
			AMP = abs(vars["ZAMP"][i])
			result = hcat(result, [Int(i), vars["TFS"][i], AMP, real(vars["ZAMP"][i]), imag(vars["ZAMP"][i])])
		end#if
    end

	return result
	end # let




end #naff

end # module
