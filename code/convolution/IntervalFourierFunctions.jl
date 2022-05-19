using IntervalArithmetic:Interval

function verifyfft(z::Vector{T}, sign=1) where T
    n = length(z); col = 1; array1 = true
    if n==1
        Z = map(T,z)
        return Z
    else
        isrow_ = false
    end
    log2n = Int(round(log2(n))) #check dimension
    if 2^log2n ≠ n #2の倍数でない場合はエラー表示
        error("length must be power of 2")
    end
    #bit-reversal(ビットリバース)
    f = 2^(log2n-1)
    v = [0;f]
    for k = 1:log2n-1
#         f = 0.5*f
        f = f >> 1
        v = append!(v,f.+v)
    end
    z2 = zeros(n,col)
    if isa(real(z[1]),Interval)
        z2 = map(T,z2)
    end
    #zを入れ替え
    for j = 1: n
        z2[j,:] = z[v[j]+1,:]
    end
    #Danielson-Lanczos algorithm
    Z = complex(map(Interval,z2))
    Index = reshape([1:n*col;],n,col)
    
#     nmax=2^14
#     if n <=nmax
#         r_real = Array{Float64}(undef, 2^14, 1)
#         file = open("INTLAB_CONST.FFTDATA_R_real.bin", "r")
#         read!(file, r_real)

#         close(file)
#         r_imag = Array{Float64}(undef, 2^14, 1)
#         file = open("INTLAB_CONST.FFTDATA_R_imag.bin", "r")
#         read!(file, r_imag)
#         close(file)

#         d = Array{Float64}(undef, 14, 1)
#         file = open("INTLAB_CONST.FFTDATA_D.bin", "r")
#         read!(file, d)
#         close(file)
        
# #         c =r[1:Int(nmax/n):nmax]
#         dd = d[log2n]
        
# #         Phizero = zeros(n)
# #         Phi = complex(map(interval,Phizero))
#         Phi = (r_real[1:Int(nmax/n):nmax] .± dd) +  im * (r_imag[1:Int(nmax/n):nmax].± dd)
# #         end
#         if sign==-1
#             Phi = adjoint.(Phi)      
#         end
#     else
#     Phi = exp.(im*theta) # SLOW (INTLAB uses table)   
#     theta = @interval(pi) * sign * (0:n-1)/n; # division exact because n is power of 2
#     Phi = cos.(theta) + im*sin.(theta) # SLOW?
    theta = map(Interval,sign * (0:n-1)/n); # division exact because n is power of 2
    Phi = cospi.(theta) + im*sinpi.(theta) # SLOW?
#     end

    v = [1:2:n;]
    w = [2:2:n;]
    t = Z[w,:]
    Z[w,:]  = Z[v,:] - t
    Z[v,:]  = Z[v,:] + t
    for index　in 1: (log2n-1)    
        m = 2^index
        m2 = 2*m
        vw = reshape([1:n;],m2,Int(n/m2))
        v = vw[1: m, :]
        w = vw[m+1: m2, : ]
        indexv = reshape(Index[v[:],:],m,Int(col*n/m2))
        indexw = reshape(Index[w[:],:],m,Int(col*n/m2))
        Phi1 = repeat(Phi[1:Int(n/m):end],outer=[1,Int(col*n/m2)])
        t = Phi1 .*  Z[indexw]
        Z[indexw] = Z[indexv] - t 
        Z[indexv] = Z[indexv] + t
    end
    reverse(Z[2:end,:],dims=2)
     if sign==-1
        Z = Z/n
    end
    if isrow_
        Z = transpose(Z)　#転置
    end
    if array1
        Z = Z[:,1]
    end
    return Z
end

function powerconvfourier(a::Vector{Complex{Interval{T}}},p) where T
    M = Int((length(a)+1)/2) # length(a) = 2M-1
    N = (p-1)*M
    ia = map(Interval, a)

    length_ia = 2*p*M-1
    length_ia_ext = nextpow(2,length_ia)# 2pM-2+2L
    
    L = Int((length_ia_ext - length_ia + 1)/2)
    
    # step.1 : padding (p-1)M + L zeros for each sides
    ia_ext = map(Complex{Interval},zeros(length_ia_ext))
    ia_ext[L+N+1:end-L-N+1] = ia  #\tilda{a}

    # step.2 : inverse fft
    ib_ext = verifyfft(ifftshift(ia_ext), -1) #sign = -1 : ifft
    
    # step.3 : power p elementwisely
    ib_extᵖ = ib_ext.^p
    
    # step.4 : fft with rescaling
    ic_extᵖ = fftshift(verifyfft(ib_extᵖ, 1)) * length_ia_ext^(p-1)  #sign = 1 : fft
    
#     return ic_extᵖ,ic_extᵖ
    return ic_extᵖ[L+N+1:end-N-L+1], ic_extᵖ[L+p:end-(L+p-2)] # return (truncated, full) version
end
