
# install.packages("Matrix")
library(Matrix)

## for estimating K (i.e. rank of E(A)):
install.packages("gdim")
library(gdim)

## for vintage sparse pca (the varimax thing):
install.packages("vsp")
library(vsp)

install.packages("irlba") # library(RSpectra) is also great! both get you sparse svd. RSpectra::svds
library(irlba)

## tidyverse is for plotting and stuff.
## This one might take awhile...
# install.packages("tidyverse")
library(tidyverse)



# this line loads the data:
load(url("https://github.com/karlrohe/dogs-say-meow/raw/master/A_journal.RData"))

# these next two lines make the graph symmetric/undirected and binary:
A = A_journal + t(A_journal)
A@x[] = 1 


#######################################################
##### ESTIMATE K WITH CROSS-VALIDATED EIGENVALUES #####
#######################################################
# https://arxiv.org/abs/2108.03336
# set.seed(1)
# maybe takes 30 seconds:
eigen_cross_validated = gdim::eigcv(A, 
                                    k_max = 100,
                                    num_bootstraps = 1,
                                    test_portion = .1,
                                    laplacian = TRUE,
                                    regularize = TRUE)

# these are z-scores (whoa):
plot(eigen_cross_validated)

# we are going to use the "regularized" graph Laplacian:
L_reg = gdim:::glaplacian(A)


##################################################
##### RADIAL STREAKS IN PCs/SINGULAR VECTORS #####
##################################################
# https://academic.oup.com/jrsssb/article/85/4/1037/7221295

# maybe takes 30 seconds:
s = irlba(L_reg,100)

# plot the top 10 singular vectors of regularized Laplacian:
u_top = s$u[,1:10]
n = nrow(u)
plot_high_leverage_sample = sample(n,1000,prob = rowSums(u_top^2))
pairs(u_top[plot_high_leverage_sample,], pch = ".")

# now let's plot 91:100:
u_bottom = s$u[,91:100]
plot_high_leverage_sample = sample(n,1000,prob = rowSums(u_bottom^2))
pairs(u_bottom[plot_high_leverage_sample,], pch = ".")


########################################################
##### RADIAL STREAKS AFTER VARIMAX ALIGN WITH AXES #####
########################################################
# https://academic.oup.com/jrsssb/article/85/4/1037/7221295

# maybe takes 30 seconds:
varimax_rotation = varimax(s$u, normalize = FALSE)$rotmat
z=s$u%*%varimax_rotation

make_z_skew_positive = sign(colSums(z^3))
z = z %*% diag(make_z_skew_positive)


z_top = z[,1:10]
plot_high_leverage_sample = sample(n,1000,prob = rowSums(z_top^2))
pairs(z_top[plot_high_leverage_sample,], pch = ".")

z_bottom = z[,91:100]
plot_high_leverage_sample = sample(n,1000,prob = rowSums(z_bottom^2))
pairs(z_bottom[plot_high_leverage_sample,], pch = ".")


################################################
##### NAME THE FACTORS / CLUSTERS / BLOCKS #####
################################################
# https://academic.oup.com/jrsssb/article/85/4/1037/7221295


uniqueJournals = rownames(A)

#bff:
source("https://raw.githubusercontent.com/karlrohe/dogs-say-meow/master/bff_with_unigrams.R")
factor_names = bff_with_unigrams(z, uniqueJournals,num_best = 5)
factor_names
# this is a lot to interpret!
# how does it all fit together????
View(factor_names) 




#################################
##### COMPUTE THE HIERARCHY #####
#################################
# https://arxiv.org/abs/2309.01301

# for this, we are going to use the asymmetric version
#   ... how a journal is cited (column information)
#   is much more informative 
#   this takes care of "spam journals" that do a lot of citing, but don't get cited much


L_reg_asymmetric = gdim:::glaplacian(A_journal)
sa= irlba(L_reg_asymmetric,100)
varimax_rotation = varimax(sa$v, normalize = FALSE)$rotmat
y=sa$v%*%varimax_rotation

make_y_skew_positive = sign(colSums(y^3))
y = y %*% diag(make_y_skew_positive)
factor_names = bff_with_unigrams(y, uniqueJournals,num_best = 5)





y_no_negs = y
y_no_negs[y_no_negs<0]=0
hist(y_no_negs[y_no_negs!=0])

B = (t(y_no_negs) %*%t(L_reg_asymmetric) ) %*% (L_reg_asymmetric %*% y_no_negs)
B = as.matrix(B)
B[B<0]=0
B = B + quantile(B[B>0],.01)
image(B)
hist(B)
S = sqrt(diag(B))
B_norm = diag(1/S)%*%B%*%diag(1/S)
B_norm[B_norm>1]=1
image(B_norm)


dist_mat= -log(B_norm)
colnames(dist_mat) = factor_names$word1
rownames(dist_mat) = factor_names$word1
library(ape)
tree_hat = ape::nj(dist_mat)
plot(tree_hat, "u", cex = .6)

