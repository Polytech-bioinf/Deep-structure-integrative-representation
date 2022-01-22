from snf import datasets
import snf

digits = datasets.load_digits()
digits.keys()

affinity_networks = snf.make_affinity(digits.data, metric='euclidean', K=20, mu=0.5)
fused_network = snf.snf(affinity_networks, K=20)
print(fused_network)