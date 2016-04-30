### POSE response

Thank you for the detailed reviews and for the opportunity to revise my paper.  I believe that I have been able to address all concerns  and that the paper is  substantially improved.  I worked hard to add the requested results and information without adding significantly to the article length, and the main paper has grown by only two pages.  There is an extensive supplemental appendix.

The dominant sentiment from the AE and the referees was that my simulation study needs to include a wider variety of error metrics and, in particular, more data generating processes.  The union of requested additional simulation configurations included smaller samples, strictly sparse models, and continuous covariates.  Each of these requests seemed to me to have merit, and since this is the Journal of _Computational_ and Graphical Statistics I figured that I should simply run everything.  

Page 15 describes all of the configurations, leading to a total of 288 different models.  I simulate and estimate 1000 times for each.  In addition to prediction error, I now provide results for estimation error on the true coefficients, the estimated model dimension, and false discovery and sensitivity rates. Incorporating results for each algorithm, and for a variety of IC and CV selection rules, leads to a set of 128 tables.  These are included in the supplemental material.  In the main paper, I have distilled results into two tables by averaging results across covariate designs but maintaining disctinction between sparsity levels, sample size, and signal strength.  Table 1 shows predictive performance, and Table 2 summarizes estimation for our sparse models (where covariate recovery is a well-posed task).  This is followed by a series of remarks that summarize the key results.

I hope that this very extensive set of simulations satisfies the request for further experimentation.   The results are broadly the same as before, with the GL procedures providing a fast alternative to the more expensive exact concave-penalized estimators.  The GL results are similar to those of our gold-standard MCP algorithm; in particular, GL with cross validation is always within 1% of MCP for predictive RMSE.  However, the more extensive simulations give a better picture of the limitations of all of the algorithms.  See pages 20-22 for detail, but some brief highlights are: 

* AICc and CV perform similarly when the sample size is large relative to the oracle support, but you need to choose carefully when the sample is too small; 

* all of the methods perform fairly well in estimation except for the marginal adaptive lasso, which does terribly;

* GL reduces false discovery relative to the lasso, but MCP manages to do so with lower drop in sensitivity

Detailed response to the AE and referee reports follows below.  Thank you again for your consideration.

##### AE report

_1. Simulation experiment is quite weak. It is a bit unusual to consider a non-sparse setup for
high-dimensional regression. It also makes little sense to compare the proposed estimator
with (18). Why not simply comparing with the true oracle? I also agree with Referee I that it
is a bit strange to choose binary regressors. By considering an AIC-type criterion, I guess
the proposed method may be useful in prediction but less in variable selection. I think it will
be more interesting to see the simulation results after selection of \\gamma._

As described above, this has all been addressed in the new simulations.  I have also added the GL-select procedure, which does as you suggest by selecting across \\gamma values.  Thank you for this idea: it works quite nicely and gives additional support for the procedure.

_2. Some more discussion on Theorem 3.1 is needed. While it compares the proposed
estimator with some L0 penalized estimator, it is not immediately clear to the reader why (6)
is a nice property without providing some discussion on the L0 method. It also does not
provide a clear guideline on how to choose the tuning parameters. The author could also
discuss whether a similar result can or cannot be achieved by Lasso._

I now preface Section 3.1 (starting page 8) with a three-paragraph discussion on why I've chosen to compare to an L0 comparator.  Also, I now reference a large literature of existing results on estimation error (of the sort that Ref I desires), provide justification for the L0 oracle as a desirable comparator, and reference later material on tuning parameter selection.

_The author seems to assume that the dispersion parameter (such as \\sigma^2 in normal
regression) is known. How if it is unknown? It is also more practical to consider \\sigma^2 to
be unknown in the simulation experiment._

I have not made such an assumption, and we always consider \\sigma unknown to the statistician.  In simulations I treat it as known for the oracle only, and in discussion of Theorem 3.1 I use it only to define an unattainable Cp oracle.  I have gone through any mention of \\sigma in the paper to try and re-phrase statements to avoid this confusion.

_The author should elaborate more on how the gamma Lasso penalty can approximate the
L0 penalty to avoid possible confusion as raised by Referee II._

This is done and I provide more description in my response to Ref II.

_The notations for the column vector and the row vector of X are quite similar. They may be
modified to avoid confusion._

Fixed; I now use \\chi instead of x for the columns.

##### Ref I

_1. The author assumes that x is replaced by x /sd(x). I understand he wanted to remove
the effect of units of x. However, from my point of view, this is not the best assumption to
reach such goal. A more common way is to assume that each x is scaled by its population
standard deviation not sample standard deviation. They are different since x/sd(x) and
x/sd(x) will be correlated with each other._

This is the only criticism that I have not addressed in actual revision.  In my theoretical setup everything is restricted to the finite sample at hand and there is no role for a population standard deviation (anyways, it would not make any difference to the results).  And in practice, of course, the population standard deviation is never available; using it in the simulation experiment would be an unrealistic advantage.

_2. The initial regularization parameter is defined as ... this is a theoretical 
definition but how do you compute the exact value in practice? I guess you assume
all weights are equal to one in the first step so it corresponds to the zero estimation in lasso
problem?_

Yes, this is easy to calculate and I now write it explicitly on page 5.  Not that it is not an assumption that all weights are one to start; this is in the algorithm definition.

_3. I assume \\hat S is the estimated support corresponding to \\lambda\_t, right? You should give the definition
there. My understanding is that the coefficient-specific weight is a decreasing function of the
coefficient estimated in the previous path step. So why not \\lambda\_{t-1}?_

Thanks, I located a typo in the algorithm that was leading to confusion.  Everything should be clear now, with \\hat S as you describe.  The weights are indeed a function of \\S\_{t-1}.

_5. It was claimed that your method can avoid biasness problem since you have different weight
for each \\beta. However, it is unclear to me that where the improvement is in Theorem 3.1. Does
your method have the problem of overfitting and shrinkage like Lasso?_

The improvement is due to small weights on true-large coefficients and large weights (near one) on the true zero coefficients.  I now remark on how the same theorem applies to the lasso, and hopefully this adds some intuition.

_ 6. Since your one-step regularization method penalizes on the model complexity, why didn’t you
come up with any simulation examples with sparse structure? And you should consider the
false positives and false negatives of your method in such setting.

7. I have the impression that your method can be very similar to adaptive lasso and I know you
discussed it in your paper. Why didn’t you also include adaptive lasso and the one step estima-
tor (Zou and Li, 2008) for comparison?

8. The sample size of you simulation example is 1000, the same as the dimensionality. This is a
very easy setting. It is much more common to consider n = 100 and p = 1000 or even more
challenging settings.

9. I am not sure why you chose binary covariates in your simulation. I suggest you consider
some cases when the covariates are continuous. In addition, the correlation structure may be
diminished since x j ’s are correlated to each other through z j ’s._

All of these concerns have been addressed in the new simulation experiment.  Note that I have always included an adaptive lasso comparator.



