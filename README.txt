This .zip file contains the R functions used for the paper
"A bivariate finite mixture growth model with selection" by
- D. Aristei (University of Perugia, IT)
- S. Bacci (University of Florence, IT)
- F.Bartolucci (University of Perugia, IT)
- S.Pandolfi (University of Perugia, IT)  

est_biv_LT.R 			estimate the finite mixture bivariate latent trajectory 
				model on the basis of an EM algorithm with NR acceleration

lk_sel.R			compute the incomplete data log-likelihood
 	
lk_sel_comp.R 			compute the complete data log-likelihood

sc_sel.R			compute the score function on the basis of the incomplete
 				data log-likelihood

sc_sel_com.R			compute the score function on the basis of the complete
				data log-likelihood

ExampleData.RData 		workspace file containing a simulated dataset


example_LT_sim.R		example file that loads the workspace ExampleData.RData					and performs model estimation

ResultsSimData.RData		workspace file containing the estimation results

