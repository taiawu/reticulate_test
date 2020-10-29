import numpy as np

def analyze_data(rdf):
  data_array = [rdf[k] for k in rdf.keys() if k != 'well']
  data_array = np.vstack(data_array)
  print(data_array.T.shape)
	  
