# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.3 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
library(cellWise)
Rcpp::sourceCpp("CppFuncs.cpp")
source("Functions.R")
source("NBA/cellMapCol.R")
load("NBA/NBAms.RData")

# https://www.nba.com/stats/players/bio?Season=2024-25&SeasonType=Regular+Season
bio = read.csv("NBA/bio.csv") |> janitor::clean_names() |> na.omit()
bio$height_cm = vapply(bio$height, to_cm, numeric(1))
bio = bio |> dplyr::select(-c(team,college, country, draft_year, draft_round, draft_number,reb,ast,height,oreb,dreb))
bio = rename(bio, c("ast" = "ast_2"))

# https://www.nba.com/stats/players/advanced?Season=2024-25&SeasonType=Regular+Season
gen = read.csv("NBA/general_adv.csv") |> janitor::clean_names() |> na.omit()
gen = gen |> dplyr::select(-c(team,age, gp, w, l,ts, offrtg, defrtg,netrtg, ast_to, ast, ast_ratio, oreb, dreb, ts, usg, pie, pace)) # pie includes ast in the formula.
d = inner_join(gen, bio, by = "player")

upper_q = quantile(as.numeric(d$height_cm))[4]

dq = d[d$height_cm >= upper_q & d$gp >= 82/2,] |> dplyr::select(-c(player,gp))
dd = d[d$height_cm >= upper_q & d$gp >= 82/2,]
set.seed(1)
A = expand.grid(rep(list(c(0, 1)), ncol(dq)-1)) |> as.matrix()
y = dq$ast
X = dq  |> dplyr::select(-c(ast)) |> scale() |> as.matrix() 
rownames(X) = NULL



set.seed(1)
ast = y
DDCdatafr = cellWise::DDC(cbind(X,ast))
fig4 = cellMapCol(t(DDCdatafr$stdResid),
           base_size = 15,
           sizetitles = 2) + xlab("Observations") + ylab("Variables") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 8))) + 
  ggtitle(element_blank())
fig4

NBAOutlier = cbind(X,y)
colnames(NBAOutlier) = c(colnames(X), "ast")
# write.csv(NBAOutlier, "NBAOutlier.csv")



start = Sys.time()
fammO = msFAMMcpp(X = X,
                  y = y,
                  A = A,
                  m = nrow(X)/2,
                  strata = 8,
                  B = 100,
                  params = "SM",
                  MMweights = T)
TimeFamm = Sys.time() - start

start = Sys.time()
frbO = msFRB.C(X = X,
               y = y,
               A = A,
               m = round(nrow(X)/2),
               B = 100,
               params = "SM",
               MMweights = T)
TimeFRB = Sys.time() - start

start = Sys.time()
mmO = msMMcpp(X = X,
              y = y,
              A = A,
              m = round(nrow(X)/2),
              strata = 8,
              B = 100,
              params = "SM",
              MMweights = T)
TimeMM = Sys.time() - start


famm0 = fammO |> dplyr::select(-c(Criterion, Time, Method))
fammModel = A[which(famm0[1,] == 1),] == 1

mm0 = mmO |> dplyr::select(-c(Criterion, Time, Method))
mmModel = A[which(mm0[1,] == 1),] == 1


frb0 = frbO |> dplyr::select(-c(Criterion, Time, Method))
frbModel = A[which(frb0[1,] == 1),] == 1


start = Sys.time()
a = computeAICcpp(A, X, y, 100)
TimeAIC = Sys.time() - start
AICModel = A[which.min(a),] == T

# Models
colnames(X)[AICModel]
colnames(X)[frbModel]
colnames(X)[fammModel]
colnames(X)[mmModel]

set.seed(1)
numFolds = 5
repeats = 10
fammErrMAE = mmErrMAE = frbErrMAE = AICErrMAE =  c()

for (rep in 1:repeats){
  folds = caret::createFolds(y, k = numFolds)
  for (i in 1:numFolds){
    Xtest = X[folds[[i]],]
    ytest = y[folds[[i]]]
    Xtrain = X[-folds[[i]],]
    ytrain = y[-folds[[i]]]
    Xytrain = cbind(ytrain,Xtrain) |> data.frame()
    
    
    fammLM  = lm(reformulate(X[,fammModel, drop = F] |> colnames(), response = "ytrain"), data = Xytrain)
    mmLM    = lm(reformulate(X[,mmModel, drop = F] |> colnames(), response = "ytrain"), data = Xytrain)
    frbLM   = lm(reformulate(X[,frbModel, drop = F] |> colnames(), response = "ytrain"), data = Xytrain)
    AICLM   = lm(reformulate(X[,AICModel, drop = F] |> colnames(), response = "ytrain"), data = Xytrain)
    
    fammPred  = predict(fammLM, data.frame(Xtest[,fammModel, drop = F]))
    mmPred   = predict(mmLM, data.frame(Xtest[,mmModel, drop = F]))
    frbPred   = predict(frbLM, data.frame(Xtest[,frbModel, drop = F]))
    AICPred   = predict(AICLM, data.frame(Xtest[,AICModel, drop = F]))
    
    fammMAE  = Metrics::mae(ytest, fammPred)
    mmMAE    = Metrics::mae(ytest, mmPred)
    frbMAE   = Metrics::mae(ytest, frbPred)
    AICMAE   = Metrics::mae(ytest, AICPred)
    
    fammErrMAE = c(fammErrMAE, fammMAE)
    mmErrMAE = c(mmErrMAE, mmMAE)
    frbErrMAE = c(frbErrMAE, frbMAE)
    AICErrMAE = c(AICErrMAE, AICMAE)
  }
  
}



mean(frbErrMAE) |> round(x = _, digits = 2)
mean(mmErrMAE) |> round(x = _, digits = 2)
mean(fammErrMAE) |> round(x = _, digits = 2)
mean(AICErrMAE) |> round(x = _, digits = 2)

sd(frbErrMAE) |> round(x = _, digits = 2)
sd(mmErrMAE) |> round(x = _, digits = 2)
sd(fammErrMAE) |> round(x = _, digits = 2)
sd(AICErrMAE) |> round(x = _, digits = 2)

save(fammErrMAE, mmErrMAE, frbErrMAE, AICErrMAE,
     TimeFamm, TimeMM, TimeFRB, TimeAIC,
     fammModel, mmModel, frbModel, AICModel,
     file = "NBAms.RData")
