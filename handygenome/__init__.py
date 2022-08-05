"""
def import_helper(*args):
	import importlib
	mod = importlib.import_module('.'.join((__name__,) + args))

	return mod
"""
