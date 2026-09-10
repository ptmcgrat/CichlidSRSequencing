import datetime, pdb, time

class Timer:
	def __init__(self):
		pass
	def test(self):
		self.start('Testing')
		time.sleep(10000)
		self.stop()
	def start(self, codeblock):
		self.start_t = datetime.datetime.now()
		print(codeblock + ': Starting....', end = '', flush = True)
	def stop(self):
		self.stop_t = datetime.datetime.now()
		delta_t = (self.stop_t - self.start_t).total_seconds()/60.0
		print('Completed in ' + '{:.1f}'.format(delta_t) + ' minutes')
