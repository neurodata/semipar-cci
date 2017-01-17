###1. Current model (joint embedding)

For group $j$ and subject $i$:

$$A_{ji} \sim Bern(logit( F C_{ji} F^T))$$
$$C_{ji} \sim N( C_j, I\sigma_1^2)$$
$$C_j \sim N(C, I\sigma_2^2)$$

Then the graphon with batch effect removed is:

$$logit( F (C_{ji}-(C_j-C)) F^T)$$

Results:

$(C_j-C)$ is really small so there weren't enough batch effects captured.



###2. Proposed model

For group j and subject i:

$$A_{ji} \sim {Bern}({logit}( F_j C_{ji} F^T_j))$$
$$vec(F_j) \sim N(vec(F), I\sigma^2)$$

Then the graphon with batch effect removed is:

$$logit( F C_{ji} F^T)$$

Good properties about this model:

1. Remove vertex-wise batch effect with F, instead of on loading L.
2. Shrinakge of the error estimates on F_j, preventing overfitting with batch-wise model.
3. Closed-form solution for optimization: Polya-Gamma EM algorithm on logit, closed-form Normal on $F_j$.