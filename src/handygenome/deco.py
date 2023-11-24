import functools
import time
import inspect
import itertools
import signal

import numpy as np

import handygenome.tools as tools
import handygenome.logutils as logutils


def make_errmsg(deconame, funcname):
    return (
        f'The names of parameters given '
        f'to the decorator generator "{deconame}"'
        f'is not included in the parameter names of '
        f'the function "{funcname}".'
    )


# timeout (https://daeguowl.tistory.com/139)
def get_deco_timeout(seconds, error_message=''):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)
        
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            #signal.alarm(seconds)
            signal.setitimer(signal.ITIMER_REAL, seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return functools.wraps(func)(wrapper)

    return decorator


def deco_timer(func):
    """Print the runtime of the decorated function"""

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time
        print(run_time)

        return value

    return wrapper_timer


def get_deco_num_set(names, n):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(make_errmsg('get_deco_num_set', func.__name__))

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            n_set = sum((name in ba.arguments) for name in names)
            if n_set != n:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    #f'number of parameters set from arguments '
                    f'number of parameters being set, '
                    f'among {tuple(names)}, must be {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_num_notNone(names, n):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                make_errmsg('get_deco_num_notNone', func.__name__)
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            n_notNone = 0
            for name in names:
                set_val = ba.arguments[name]
                if set_val is not None:
                    n_notNone += 1

            if n_notNone != n:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    f'number of parameters, among {tuple(names)}, '
                    f'being set as a value different from the default, '
                    f'must be {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_num_set_differently(names, n, how='equal'):
    assert how in ('equal', 'gt', 'lt')
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                make_errmsg('get_deco_num_set_differently', func.__name__)
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            n_diff = 0
            for name in names:
                default_val = sig.parameters[name].default
                set_val = ba.arguments[name]
                if default_val is None:
                    if set_val is not None:
                        n_diff += 1
                else:
                    if set_val != default_val:
                        n_diff += 1

            if how == 'equal':
                cond = (n_diff == n)
                word = 'equal to'
            elif how == 'gt':
                cond = (n_diff > n)
                word = 'greater than'
            elif how == 'lt':
                cond = (n_diff < n)
                word = 'less than'

            if not cond:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    f'number of parameters, among {tuple(names)}, '
                    f'being set as a value different from the default, '
                    f'must be {word} {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_arg_choices(mapping):
    """Args:
        mapping: {'argname': (valid_value1, valid_value2, ...), ...}
    """

    def decorator(func):
        sig = inspect.signature(func)
        if not set(mapping.keys()).issubset(sig.parameters.keys()):
            raise Exception(
                make_errmsg('get_deco_arg_choices', func.__name__)
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            for key, val in mapping.items():
                if ba.arguments[key] not in val:
                    raise ValueError(
                        f'For the function "{func.__name__}", '
                        f'the parameter "{key}" must be one of these values: '
                        f'{tuple(val)}.'
                    )

            return func(*args, **kwargs)

        return wrapper

    return decorator

# other decorators
def get_deco_timestamp(msg, logger):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger.info(f'BEGINNING {msg}')
            result = func(*args, **kwargs)
            logger.info(f'FINISHED {msg}')
            return result
        return wrapper

    return decorator


def scalarize_result(result):
    if isinstance(result, (list, tuple)):
        result = tuple(
            x[0] if x.shape == (1,) else x
            for x in result
        )
    else:
        if result.shape == (1,):
            result = result[0]


def vectorize(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        keys = list(kwargs.keys())
        arglist = list(itertools.chain(args, kwargs.values()))
        arglist = np.broadcast_arrays(arglist)

        new_args = arglist[:len(args)]
        new_kwargs = dict(zip(keys, arglist[len(args):]))

        result = func(*new_args, **new_kwargs)
        if isinstance(result, (list, tuple)):
            result = tuple(
                x[0] if x.shape == (1,) else x
                for x in result
            )
        else:
            if result.shape == (1,):
                result = result[0]

        return result

    return wrapper


def args_into_array(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        keys = list(kwargs.keys())
        arglist = [
            np.atleast_1d(x) for x in 
            itertools.chain(args, kwargs.values())
        ]

        new_args = arglist[:len(args)]
        new_kwargs = dict(zip(keys, arglist[len(args):]))

        result = np.squeeze(func(*new_args, **new_kwargs))

        return result

    return wrapper


def get_deco_atleast1d(names, keep_none=True):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_atleast1d" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            for key in names:
                old = ba.arguments[key]
                if old is None:
                    if keep_none:
                        new = old
                    else:
                        new = np.atleast_1d(old)
                elif np.isscalar(old):
                    new = np.atleast_1d(old)
                else:
                    new = np.atleast_1d(tuple(old))
                ba.arguments[key] = new

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator


def get_deco_squeeze_atleast1d(names):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_squeeze_atleast1d" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            for key in names:
                ba.arguments[key] = np.atleast_1d(
                    np.squeeze(
                        ba.arguments[key]
                    )
                )

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator


def get_deco_broadcast(names):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                make_errmsg('get_deco_broadcast', func.__name__)
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            bc_args = np.broadcast_arrays(*[ba.arguments[key] for key in names])
            for key, newarg in zip(names, bc_args):
                ba.arguments[key] = newarg

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator


def get_deco_asarray(names):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                make_errmsg('get_deco_asarray', func.__name__)
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            for key in names:
                ba.arguments[key] = np.asarray(ba.arguments[key])

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator


def get_deco_nproc_limit(max_nproc):
    def decorator(func):
        sig = inspect.signature(func)
        if 'nproc' not in sig.parameters.keys():
            raise Exception(f'"nproc" must be included in the decorated function parameters')
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            if ba.arguments['nproc'] > max_nproc:
                logutils.log(f'Maximum allowed "nproc" is {max_nproc}. "nproc" value of {max_nproc} will be used.', level='warning')
                ba.arguments['nproc'] = max_nproc

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator



