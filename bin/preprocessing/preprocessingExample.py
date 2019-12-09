
from os.path import *
from Preprocessing import Preprocessing


processor = Preprocessing(join(dirname(__file__),"exampleFiles"))
processor.run()