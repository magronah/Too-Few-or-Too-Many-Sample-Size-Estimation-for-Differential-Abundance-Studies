# Too-Few-or-Too-Many-Sample-Size-Estimation-for-Differential-Abundance-Studies

## Abstract
Determining an appropriate sample size for a study is a crucial step in planning 
scientific research.
Appropriate sample size planning avoids both inadequate and inflated sample sizes.  
Inflated sample sizes wastes resources, time and effort of human subjects, and 
lives of experimental animals. Inadequate sample sizes, a much more common problem,
wastes even more resources through the inability to detect biologically 
meaningful differences and encourages questionable research practices 
like $p$-hacking.  
Microbiome studies are particularly challenged by small sample sizes, particularly
in studies of human subjects or expensive animal models. 
In practice, the statistical power of taxa within a differential abundance study 
is influenced by the effect size (typically quantified as fold change), 
mean abundance of individual taxa, and the number of samples. 
We present a novel approach for sample size calculation for differential 
abundance studies as a function of effect size, mean abundance and statistical
power of taxa. 
Our method is implemented in the [`power.nb`](https://github.com/magronah/power.nb) 
`R` package, available at 
available on GitHub at 
[`(https://michaelagronah.com/power.nb/articles/stub.html)`](https://michaelagronah.com/power.nb/articles/stub.html).
We applied our model for sample size calculation using estimates
of mean abundance and fold change of taxa obtained from thirty 
real-world microbiome datasets. Our results showed that differential 
abundance microbiome studies require larger sample sizes than are currently 
prevalent in the literature to achieve adequate statistical power. 
Our framework will help researchers make informed decisions about 
appropriate sample sizes.
