import time
import cProfile
import pstats
from functools import wraps
from pathlib import Path
from datetime import datetime

def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()

        print(f"{func.__name__}: {end-start:.6f}s")
        return result

    return wrapper


def profiler(logfile=None, sort_by="cumtime"):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            prof = cProfile.Profile()

            prof.enable()
            result = func(*args, **kwargs)
            prof.disable()

            if logfile:
                with open(logfile, "w") as f:
                    stats = pstats.Stats(prof, stream=f)
                    stats.sort_stats(sort_by)
                    stats.print_stats()
            else:
                stats = pstats.Stats(prof)
                stats.sort_stats(sort_by)
                stats.print_stats()

            return result

        return wrapper

    return decorator


def profile(
    logfile=True,
    save_prof=True,
    sort_by="cumtime",
    top_n=20,
    output_dir="profiles",
):
    """
    Decorator for timing and profiling functions.

    Example:
        @profile()
        def foo():
            pass
    """

    def decorator(func):

        @wraps(func)
        def wrapper(*args, **kwargs):

            Path(output_dir).mkdir(parents=True, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            log_path = (
                Path(output_dir)
                / f"{func.__name__}_{timestamp}.log"
            )

            prof_path = (
                Path(output_dir)
                / f"{func.__name__}_{timestamp}.prof"
            )

            profiler = cProfile.Profile()

            start = time.perf_counter()

            profiler.enable()
            result = func(*args, **kwargs)
            profiler.disable()

            elapsed = time.perf_counter() - start

            if save_prof:
                profiler.dump_stats(str(prof_path))

            if logfile:
                with open(log_path, "w") as f:

                    f.write(
                        f"Function : {func.__name__}\n"
                    )
                    f.write(
                        f"Runtime  : {elapsed:.6f} s\n\n"
                    )

                    stats = pstats.Stats(
                        profiler,
                        stream=f
                    )

                    stats.strip_dirs()
                    stats.sort_stats(sort_by)
                    stats.print_stats(top_n)

            print(
                f"[PROFILE] "
                f"{func.__name__}: "
                f"{elapsed:.6f}s"
            )

            return result

        return wrapper

    return decorator