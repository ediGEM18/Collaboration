import pickle

def read(filename):
	file = open(filename, 'rb')
	return pickle.load(file)

def write(filename, to_pickle):
	file = open(filename, 'wb')
	pickle.dump(to_pickle, file)
	file.close()
