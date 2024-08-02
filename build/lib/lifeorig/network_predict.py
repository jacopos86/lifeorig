import warnings
from lifeorig.logging_module import log
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from tensorflow import keras
#
#  This class defines a predictor
#  network for the fitness
#  1) multilayer perceptron
#  2) keras DNN -> pytorch model
warnings.filterwarnings("ignore")
#
#  driver function to create
#  the network model
def generate_predictor(NN_model, size):
    if NN_model == "MLP":
        return MLP_predictor(size)
    elif NN_model == "DNN":
        return DNN_predictor(size)
    else:
        log.error("Wrong NN model in input")
#  build X, y input data
#  for the network
def build_Xy_data(ACF_set):
    X = []
    y = []
    for ACF in ACF_set:
        print(len(ACF.genome))
        inp = []
        for c in ACF.genome:
            inp.append(int(c))
        print(inp)
        X.append(inp)
        y.append(ACF.fitness)
    return X, y
#
# base class
class fitness_predictor:
    def __init__(self, genome_size):
        self.inp_size = genome_size
        self.model = None
        self.regr = None
#
#  MLP class
class MLP_predictor(fitness_predictor):
    def __init__(self, genome_size):
        super(MLP_predictor, self).__init__(genome_size)
    # accumulate network data
    # this is used for setting model
    def set_model(self, NN_parameters):
        n_hidden_layers = tuple(NN_parameters['n_hidden_layers'])
        activation = NN_parameters['activation']
        solver = NN_parameters['solver']
        alpha = NN_parameters['alpha']
        max_iter = NN_parameters['max_iter']
        random_state = NN_parameters['random_state']
        # define MLP model
        self.model = MLPRegressor(hidden_layer_sizes=n_hidden_layers,
                    activation=activation, solver=solver,
                    alpha=alpha, random_state=random_state,
                    max_iter=max_iter)
    # fit data -> model training
    def fit(self, NN_parameters, X, y):
        random_state = NN_parameters['random_state']
        test_size = NN_parameters['test_size']
        shuffle = NN_parameters['shuffle']
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=random_state, test_size=test_size, shuffle=shuffle)
        # build regressor
        self.regr = self.model.fit(X_train, y_train)
        return X_test, y_test
    # get score
    def get_score(self, X_test, y_test):
        log.info("\t N. LAYERS MULTILAYER PERCEPTRON MODEL : " + str(self.regr.n_layers_))
        log.info("\t MODEL SHAPE : " + str(len(self.regr.coefs_)))
        score = self.regr.score(X_test, y_test)
        return str(score)
    # predict value
    def predict(self, X):
        y = self.regr.predict(X)
        return y
#
# DL class
class DNN_predictor(fitness_predictor):
    def __init__(self, genome_size):
        super(DNN_predictor, self).__init__(genome_size)
    # accumulate network data
    # this is used for setting model
    def set_model(self, NN_parameters):
        n_hidden_layers = NN_parameters['n_hidden_layers']
        n_hidden_units = NN_parameters['n_hidden_units']
        activation_in = NN_parameters['activation_in']
        activation_hid = NN_parameters['activation_hid']
        activation_out = NN_parameters['activation_out']
        loss = NN_parameters['loss']
        optimizer = NN_parameters['optimizer']
        input_shape = self.inp_size
        # build the model
        self.model = keras.Sequential()
        self.model.add(keras.layers.Dense(units=1, activation=activation_in, input_shape=[input_shape]))
        for n in range(n_hidden_layers):
            self.model.add(keras.layers.Dense(units=n_hidden_units, activation=activation_hid))
        self.model.add(keras.layers.Dense(units=1, activation=activation_out))
        # compile model
        self.model.compile(loss=loss, optimizer=optimizer, metrics=['accuracy'])
        # display model
        info = str(self.model.summary())
        log.info(info)
    # fitting NN model
    def fit(self, NN_parameters, X, y):
        epochs = NN_parameters['epochs']
        verbose= NN_parameters['verbose']
        random_state = NN_parameters['random_state']
        test_size = NN_parameters['test_size']
        shuffle = NN_parameters['shuffle']
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=random_state, test_size=test_size, shuffle=shuffle)
        # fitting
        self.model.fit(X_train, y_train, epochs=epochs, verbose=verbose)
        print(self.model.predict(X_train), y_train)
        return X_test, y_test
    # score
    def get_score(self, X_test, y_test):
        _, test_accuracy = self.model.evaluate(X_test, y_test)
        print(self.model.predict(X_test)[:10], y_test[:10])
        return str(test_accuracy*100.)