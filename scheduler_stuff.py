import dask
dask.config.set(scheduler='threads')

from dask.distributed import Client
if __name__ == '__main__':
    file = os.getenv('MEMBERWORK') + '/gen010/my-scheduler.json'
    client = Client(scheduler_file=file)
