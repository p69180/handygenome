import json
import re
import contextlib
import time
import shutil
import subprocess

import ftplib
import urllib.request
import urllib.error
import urllib.parse

import handygenome.logutils as logutils
import handygenome.deco as deco
import handygenome.tools as tools


HTTP_HEADER_POST = {
    "Content-Type": "application/json", 
    "Accept": "application/json",
}
HTTP_HEADER_GET = {'Content-Type': 'application/json'}


class MaxRetryError(Exception):
    pass


def retry_or_raise(n_try, retry_count, exc, retry_interval, msg_prefix=''):
    """helper function"""
    if (retry_count is None) or (n_try <= retry_count):
        retry_count_string = (
            'inf' if retry_count is None else str(retry_count)
        )
        logutils.log(
            f'{msg_prefix}Retrying due to failure (Exception: {exc}); n_try={n_try}/{retry_count_string}',
            level='info',
        )
#        if isinstance(exc, TimeoutError):
#            logutils.log(
#                f'{msg_prefix}Retrying due to TimeoutError (Exception: {exc}); '
#                f'n_try={n_try}/{retry_count_string}'
#            )
#        else:
#            logutils.log(
#                f'{msg_prefix}Retrying due to non-TimeoutError (Exception: {exc}); '
#                f'n_try={n_try}/{retry_count_string}'
#            )
        time.sleep(retry_interval)
    else:
        if isinstance(exc, urllib.error.URLError):
            print(exc.read().decode('utf-8'))
        else:
            print(str(exc))
        raise MaxRetryError(f'Exceeded maximum retry count({retry_count}).') from exc


# ftp

def ftp_login(url, retry_count=10, retry_interval=1, timeout=5):
    @deco.get_deco_timeout(timeout)
    def helper():
        ftp = ftplib.FTP(url, timeout=timeout)
        ftp.login()
        return ftp

    retry_msg_pf = 'ftp login: '

    n_try = 0
    while True:
        n_try += 1
        try:
            #with contextlib.redirect_stdout('/dev/null'):
            ftp = helper()
            #ftp = ftplib.FTP(url, timeout=timeout)
            #ftp.login()
        except TimeoutError as exc:
            retry_or_raise(n_try, retry_count, exc, retry_interval, msg_prefix=retry_msg_pf)
            continue
        except OSError as exc:
            if (str(exc) == '[Errno 101] Network is unreachable'):
                retry_or_raise(n_try, retry_count, exc, retry_interval, msg_prefix=retry_msg_pf)
                continue
            else:
                raise
        else:
            break

    return ftp



def ftp_nlst(
    authority, path, 
    retry_count=10, retry_interval=1, timeout=5, verbose=False,
):
    retry_msg_pf = 'ftp nlst: '

    if verbose:
        logutils.log(
            f'Performing ftp nlst (authority={repr(authority)}, path={repr(path)})',
            level='debug',
        )

    n_try = 0
    while True:
        n_try += 1
        try:
            ftp = ftp_login(
                authority, 
                retry_count=1, 
                retry_interval=retry_interval, 
                timeout=timeout,
            )
            result = ftp.nlst(path)
        except MaxRetryError as exc:  # ftp_login failure
            retry_or_raise(n_try, retry_count, exc, retry_interval, msg_prefix=retry_msg_pf)
            continue
        except TimeoutError as exc:  # nlst failure
            retry_or_raise(n_try, retry_count, exc, retry_interval, msg_prefix=retry_msg_pf)
            continue
#        except OSError as exc:  # nlst failure when ftp object is timed-out
#            if (str(exc) == 'cannot read from timed out object'):
#                retry_or_raise(n_try, retry_count, exc, retry_interval, msg_prefix=retry_msg_pf)
#                continue
#            else:
#                raise
        else:
            break

    return result


def trim_path(path):
    path = re.sub('/+$', '', path)
    path = re.sub('/{2,}', '/', path)
    return path


def join_ftp_paths(*args):
    return trim_path('/'.join(args))


def ftp_listdir(ftp, path):
    return ftp.nlst(trim_path(path))


def ftp_listdir_old(ftp, path):
    path = trim_path(path)
    fname_list = list()
    ftp.cwd(path)
    ftp.retrlines('NLST', callback=(lambda x: fname_list.append(path + '/' + x)))
    return fname_list


# http

def http_run_urlopen(url_or_req, retry_count=10, retry_interval=1, urlopen_timeout=5):
    url_string = (
        url_or_req
        if isinstance(url_or_req, str) else
        url_or_req.full_url
    )
    logutils.log(f'Trying to open url {repr(url_string)}', level='info')

    n_try = 0
    while True:
        n_try += 1
        try:
            response = urllib.request.urlopen(url_or_req, timeout=urlopen_timeout)
        except urllib.error.URLError as exc:
            if isinstance(exc, urllib.error.HTTPError):
                if exc.code == 500:  # internal server error
                    retry_or_raise(n_try, retry_count, exc, retry_interval)
                    continue
                else:
                    print(exc.read().decode('utf-8'))
                    raise
            else: 
                if str(exc) == '<urlopen error [Errno 101] Network is unreachable>':
                    retry_or_raise(n_try, retry_count, exc, retry_interval)
                    continue
                else:
                    print(exc.read().decode('utf-8'))
                    raise
        except TimeoutError as exc:
            retry_or_raise(n_try, retry_count, exc, retry_interval)
            continue
        else:
            break

    logutils.log(f'Succeeded to open url {repr(url_string)}', level='info')

    return response


def http_get(url, params=None, headers=None, text=False, retry_count=10, retry_interval=1):
    # set params
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)
    if headers is None:
        if text:
            headers = {'Accept': 'text/plain'}
        else:
            headers = {'Accept': 'application/json'}
    # main
    req = urllib.request.Request(url, headers=headers, method='GET')
    return http_send_request(req, text, retry_count, retry_interval)


def http_post(url, data, params=None, headers=None, text=False, retry_count=10, retry_interval=1):
    # set params
    data_code = json.dumps(data).encode('ascii')
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)
    if headers is None:
        if text:
            headers = {'Accept': 'text/plain'}
        else:
            headers = {'Accept': 'application/json'}
    # main
    req = urllib.request.Request(url, data=data_code, headers=headers, method='POST')
    result = http_send_request(req, text, retry_count, retry_interval)
    return result


def http_send_request(req, text, retry_count=10, retry_interval=1, urlopen_timeout=5):
    with http_run_urlopen(
        req, 
        retry_count=retry_count, 
        retry_interval=retry_interval, 
        urlopen_timeout=urlopen_timeout,
    ) as response:
        if text:
            result = response.read().decode('utf-8')
        else:
            result = json.loads(response.read())

    return result


def download(url, path, retry_count=10, retry_interval=1, urlopen_timeout=5):
    logutils.log(f'Beginning download: url={repr(url)}, path={repr(path)}', level='info')
    while True:
        try:
            with http_run_urlopen(
                url, 
                retry_count=retry_count, 
                retry_interval=retry_interval, 
                urlopen_timeout=urlopen_timeout,
            ) as response:
                with open(path, 'wb') as outfile:
                    shutil.copyfileobj(response, outfile)
        except TimeoutError as exc:
            logutils.log(f'Retrying due to TimeoutError (Exception: {exc})', level='info')
            continue
        else:
            break
    logutils.log(f'Finished download: url={repr(url)}, path={repr(path)}', level='info')


def download_wget(url, path):
    subprocess.run(['wget', '-O', path, url])

