library(rbacon)

# run Bacon.cleanup() if necessary

Bacon('HIBO19')

# print age and 95% confidence interval at chosen depth
print(mean(Bacon.Age.d(10.5)))
print(quantile(Bacon.Age.d(10.5), probs=c(.025, .975)))