import datetime
import time
import multiprocessing

def format_time(timestamp):
    return timestamp.strftime("%Y-%m-%d %H:%M:%S")

def worker(task_number):
    start_time = datetime.datetime.now()
    print(f"[{format_time(start_time)}] Task {task_number} started.")
    
    # Simulate some work with sleep
    time.sleep(task_number)
    
    end_time = datetime.datetime.now()
    print(f"[{format_time(end_time)}] Task {task_number} finished.")

if __name__ == "__main__":
    # Number of concurrent processes
    max_processes = 3
    # Example tasks
    tasks = [5, 3, 8, 2, 6, 1, 7, 4]

    with multiprocessing.Pool(processes=max_processes) as pool:
        pool.map(worker, tasks)
