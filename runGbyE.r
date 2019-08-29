runGbyE <- function(dataset, Y = "Y", G = "G", E = "E", covars = NULL, funct = "binomial(link='logit')", iterm = "I(G * E)") {
    require("sandwich", quietly = T)
    require("mgcv", quietly = T)
    require("data.table", quietly = T)
    
    dataset <- data.table(dataset)
    
    covterm <- ifelse(!is.null(covars), paste0("+", paste(covars, collapse = "+")), "")
    dataset <- na.omit(dataset[, c(Y, G, E, covars), with = F])
    setnames(dataset, c(Y, G, E), c("Y", "G", "E"))
    
    mod1 <- eval(parse(text = paste("glm(Y~G+E+", iterm, covterm, ", data=dataset, family=", funct, ")")))

    ## Model-based P
    pval.model <- summary(mod1)$coefficients["I(G * E)", 4]
  
  	## Sandwich robust variance-based P  
    sandwich_se1 <- diag(vcovHC(mod1, type = "HC"))^0.5
    z1 <- coef(mod1)/sandwich_se1
    pval.sandwich <- as.numeric(pchisq(z1^2, 1, lower.tail = FALSE)["I(G * E)"])  #get the pvalue corresponding interaction
    
    ## GAM
    mod2 <- eval(parse(text = paste0("gam(Y ~ s(E) + G +", iterm, covterm, ", data=dataset, family=", funct, ")")))
    res.gam <- summary(mod2)
    pval.gam <- unname(res.gam$p.pv[iterm])
    
    return(data.table(G = G, `Model-based P` = pval.model, `Sandwich robust variance-based P` = pval.sandwich, `GAM P` = pval.gam))
}
