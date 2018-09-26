


import bob.learn.em
import bob.db.iris
import numpy
import matplotlib.pyplot as plt
import sklearn
#from sklearn import preprocessing
#import gmm_top_helpers as gmm_t
#from sklearn.model_selection import GroupShuffleSplit
#from sklearn import preprocessing

print('The scikit-learn version is {}.'.format(sklearn.__version__))


data_per_class = bob.db.iris.data()
setosa = numpy.column_stack(
    (data_per_class['setosa'][:, 0], data_per_class['setosa'][:, 3]))
versicolor = numpy.column_stack(
    (data_per_class['versicolor'][:, 0], data_per_class['versicolor'][:, 3]))
virginica = numpy.column_stack(
    (data_per_class['virginica'][:, 0], data_per_class['virginica'][:, 3]))

data = numpy.vstack((setosa, versicolor, virginica))

# Two clusters with a feature dimensionality of 3
machine = bob.learn.em.GMMMachine(3, 2)
trainer = bob.learn.em.ML_GMMTrainer(True, True, True)
machine.means = numpy.array([[5, 3], [4, 2], [7, 3.]])
bob.learn.em.train(trainer, machine, data, max_iterations=200,
                   convergence_threshold=1e-5)  # Train the KMeansMachine