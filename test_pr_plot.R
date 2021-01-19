library(reticulate)
#use_virtualenv("nmf")

source_python("read_pickle.py")

fn = '/Users/liuya/Downloads/ctimage/01-18/MethPerf-K562_WGBS.plot.curve.data.ytrue.ypred.Singleton.pkl'

ret <- read_pickle_file(fn)

#length(ret['DeepSignal_true'])
#length(ret['DeepSignal_pred'])

ret