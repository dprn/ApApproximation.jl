##############################
# BISPECTRAL SETS GENERATION #
##############################

def_rho = 84
N = 64
camambert_elem = 4
def_theta = N*camambert_elem

tot_lambda = 85
max_lambda = 2. #1.25 # 1.5 seems to be the Shannon limit
def_lambda = max_lambda/tot_lambda
camambert_lambda = 4

cutoff = [	(camambert_elem,1:1.:def_rho/2),
			(camambert_elem,def_rho/2+1:1.:def_rho)
		];
E = BispectralSet(N, cutoff)

cutoff = [ (camambert_lambda,def_lambda:def_lambda:tot_lambda*def_lambda/2),
    (camambert_lambda,tot_lambda*def_lambda/2+def_lambda:def_lambda:tot_lambda*def_lambda)
		];
def_theta_lambda = N*camambert_lambda
F = BispectralSet(N, cutoff);

###################
# BESSEL MATRICES #
###################

# Weights
g1 = .1
g2 = .1
g3 = 1.
g4 = 1e2
morceau = round(Int,size(F,1)/4)
D = [ g1*ones(morceau); g2*ones(morceau); g3*ones(morceau); g4*ones(size(F,1) - 3*morceau); ]
D *= 100

println("* Generating Bessel matrices")
@time J = BesselMatrix(E,F, weights = D);

####################
# AP APPROXIMATION #
####################

test_img = "lena_gray_256"
img = convert(Array{Float64,2}, testimage(test_img) |> data)

println("* Bispectral interpolation")
@time f = cartesian2bispectral(img, E)

println("* AP interpolation")
@time af = ap(f, F, J)

println("* Rotation and inverse AP")
rot = 36/120 * pi
af_rot = rotate(af, rot)
@time f_rot = iap(af_rot, E, J)

println("* Translation and inverse AP")
ρ = 100.
θ = pi/3
af_trans = translate(af, ρ, θ)
@time f_trans = iap(af_trans, E, J)

println("* Bispectral to Cartesian conversion")
def_x = 256
def_y = 256
bisp2cart(f) = pol2cart(bispectral2pol(f), def_x, def_y, clip = false)
res_rot= bisp2cart(f_rot)
res_trans = bisp2cart(f_trans)

###############
# SAVE IMAGES #
###############

save_result(f, file) = save("$file.png", convert(Image, clamp(f,0,1)))
save_result(res_rot, "rotation")
save_result(res_trans, "translation")
