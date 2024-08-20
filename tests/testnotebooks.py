__author__ = 'Robert Nikutta <robert.nikutta@noirlab.edu>, Data Lab Team <datalab@noirlab.edu'
__version__ = '20240819'

# imports

# stdlib
import glob
import json
import os
import re
import time
from getpass import getpass
from pprint import PrettyPrinter
pp = PrettyPrinter()
pprint = pp.pprint

# third party
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbclient.exceptions import DeadKernelError
from jupyter_client.kernelspec import NoSuchKernel
import pandas as pd
pd.set_option('max_colwidth', 800)

# Data Lab
from dl.authClient import login

# Compile regex to filter out ANSI escape codes
ansi_escape = re.compile(r'(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]')


def cprint(msg,color='green',bar=7,pad=' ',frame=False,newline=False):
    """Color-print msg.

    Parameters
    ----------
    msg : string
        Message to be color-printed.

    color : string
        Color to use for printing. Currently defined are 'black',
       'red', 'green' (default), 'yellow', 'blue'.

    pad : str
       String to insert between bar and msg, both left and
       right. Default: ' '

    bar : int
        Length of solid-color characters (a bar) to print before and
        after msg. Default: 7

    frame : bool
        If True, print message surrounded with a solid-color
        frame. Default: False.

    newline : bool
        If True, print one newline at the end. Default: False

    Examples
    --------

    cprint('Hello world!',color='blue',bar=7,frame=False)
    ███████ Hello world! ███████

    cprint('Hello all worlds!',color='green',bar=1,frame=True)
    █████████████████████
    █ Hello all worlds! █
    █████████████████████

    """
    
    base = 30
    colors = {'black' : 0,
              'red' : 1,
              'green' : 2,
              'yellow' : 3,
              'blue' : 4}
    
    c = colors[color]
    s = '\x1b[0;%d;%dm' % (base+c,base+10+c) # start
    t = '\x1b[0;%d;1m' % (base+c)            # text
    e = '\x1b[0m'                            # end

    l = len(msg)
    
    def _print(aux):
        print(s + aux + e)
    
    if frame is True:
        l = len(msg)
        tl = 2*bar + l + 2
        line = ' '*tl

    if frame: _print(line)
    _print(' '*bar + t + pad + msg + pad + s + ' '*bar)
    if frame: _print(line)

    if newline is True:
        print()
        

def run_notebook(nbpath,kernel='python3'):

    """Execute notebook with nbconvert, save to a new file, and parse for
errors.

    Parameters
    ----------
    nbpath : str
       Path to the .ipynb notebook file to be tested.

    kernel : str
       Name of kernel to be used. Default: python3
       Currently kernel switching is not yet implemented.

    Attribution
    -----------
    Some code in this function from
    http://www.blog.pythonlibrary.org/2018/10/16/testing-jupyter-notebooks/

    """
    
    nbname, _ = os.path.splitext(os.path.basename(nbpath))
    dirname = os.path.dirname(nbpath)
    
    # check if any code needs to be injected as the first cell in a notebook
    if os.path.isfile('prerun.py'):
        with open('prerun.py','r') as prerun:
            lines = prerun.readlines() # read on all lines
            lines = [line for line in lines if not line[0] in ('#','\n')] # remove comments and empty lines
            if len(lines) > 0:
                lines.insert(0,"# Code cell inserted from prerun.py file. Comment out in prerun.py if you don't wish to run it.\n")
                lines[-1] = lines[-1].strip() # strip line break from very last line
    else:
        lines = []
    
    with open(nbpath) as f:
        nb = nbformat.read(f, as_version=4)
        
        # inject a code cell as first cell, if there is code in prerun.py
        if len(lines) > 0:
            nb.cells.insert(0,nbformat.v4.new_code_cell(''.join(lines)))
        
    proc = ExecutePreprocessor(timeout=6000, kernel_name=kernel)
    proc.allow_errors = True
 
    proc.preprocess(nb, {'metadata': {'path': dirname}})
    output_path = os.path.join(dirname, '{}_tested.ipynb'.format(nbname))
 
    with open(output_path, mode='wt') as f:
        nbformat.write(nb, f)
 
    errors = []
    for cell in nb.cells:
        if 'outputs' in cell:
            for output in cell['outputs']:
                if output.output_type == 'error':
                    errors.append(output)
 
    return nbpath, nb, errors


def glob_with_pattern(paths,patterns=('/**/*.ipynb',)):
    
    if isinstance(paths,str):
        paths = (paths,)
        
    if isinstance(patterns,str):
        patterns = (patterns,)
            
    files = []
    
    for _ in paths:
#        if os.path.isfile(_):
#            if _.endswith('.ipynb'):
#                files.append(_)
        if os.path.isdir(_):
            for pat in patterns:
                files += glob.glob(_ + pat, recursive=True)
        else:
            print('Not a dir:', _)
            
    return set(files)


def myglob(paths,include=('/**/*.ipynb',),exclude=None,exclude_default=\
           ('/**/*_tested.ipynb','**/testnotebooks.ipynb',\
            '**/*Template.ipynb')):
    
    """Find files that match some patterns, but exclude other patterns.
   
    For paths given by `paths`, find all file path names that match any pattern
    in `include`, but exclude all file path names that match any pattern in
    `exclude` or `exclude_default`.
    
    Parameters:
    -----------    
    paths : str or tuple
        Either one (if string) or multiple (if tuple of strings) paths to scan
        for the desired file path names. `paths`, or the members of `paths`,
        can be either file name paths, or directories to search. Each directory
        will be searched recursively.
        
    include: str or tuple
        File path name pattern (if string) or tuple of patterns to match found
        file path names against. By default, all ''*ipynb' files will be found.
        
    exclude: str or tuple
        File path name pattern (if string) or tuple of patterns to be excluded
        from the final set of matched file path names. Default is None.
    
    exclude_default: str or tuple
        File path name pattern (if string) or tuple of patterns, to be excluded
        by default from the final set of matched file path names. This is in
        addition to any exclusion patterns given in `exclude`. By default,
        excluded patterns comprise: notebook files ending in '*_tested.ipynb',
        files named 'testnotebooks.ipynb', and notebook files ending in
        '*Template.ipynb'.
    
    Examples:
    ---------
    
    Find all ''*.ipynb' files in directories ''../a/b/dir1/'' and ./XYZ/dir2/'
    but exclude those ending in '_foobar.ipynb':

    ```python
    paths = ('../a/b/dir1/', './XYZ/dir2/')
    exclude = '**/*_foobar.ipynb'
    files = myglob(paths=paths,exclude=exclude)
    ```
    
    """

    if isinstance(paths,str):
        paths = (paths,)
        
    files = set([_ for _ in paths if os.path.isfile(_)])
    dirs = [_ for _ in paths if os.path.isdir(_)]
    
    # set of all files to be excluded
    files_exclude = set()
    if len(dirs) > 0:
        # set of all matched files
        files = files.union( set(glob_with_pattern(paths,patterns=include)) )
        
        if exclude is not None:
            aux = set(glob_with_pattern(dirs,patterns=exclude))
            files_exclude = files_exclude.union(aux)

        if exclude_default is not None:
            aux = set(glob_with_pattern(dirs,patterns=exclude_default))
            files_exclude = files_exclude.union(aux)
        
    # all-excluded
    return sorted(files - files_exclude)  # set operation



def get_kernel_name(nbpath,default='python3'):
    with open(nbpath,'r') as f:
        aux = json.load(f)

    try:
        kernel = aux['metadata']['kernelspec']['name']
    except KeyError:
        kernel = default

    return kernel

def get_nbs(paths,include=('/**/*.ipynb',),exclude=('/**/*_tested.ipynb',)):

    """Construct a list of all *.ipynb files to run.

    Run all notebooks found in `paths`, excluding those that match name patterns given by `exclude`.

    Parameters:
    -----------

    paths: str or tuple
        See docstring of :func:`myglob`.

    include: str or tuple
        See docstring of :func:`myglob`.

    exclude: str or tuple
        See docstring of :func:`myglob`.

    exclude_default: str or tuple
        See docstring of :func:`myglob`.


    Returns:
    --------

    nbs : list
        List of valid notebook files to test.
    """

    nbs = myglob(paths,include=include,exclude=exclude)
    
    print()
    cprint('Will test these notebooks:',bar=0,color='yellow',frame=False,pad='')
    pprint(nbs)
    cprint('Number of notebooks: ' + str(len(nbs)),bar=0,color='yellow',frame=False,pad='',newline=True)
    
    return nbs


def run(nbs,plain=False):

    """Run *.ipynb files given in the `nbs` list via nbconvert and report PASS/FAIL test matrix.

    Also record runtime for each notebook file.

    Parameters:
    -----------

    nbs: list
        List of notebook files (full paths) to run in the test.

    plain: bool
        If False (the default, for usage in Jupyter and color-aware
        terminals), report errors with ANSI escape codes for color. If
        True, filter out these escape codes (i.e. print in plain text)

    """
    
    tests = []  # will record test results per notebook file
        
    for j, nbfile in enumerate(nbs):
        start = time.time()

        print('================================================')
        print('TESTING NOTEBOOK %d/%d: %s' % (j+1,len(nbs),nbfile))

        kernel = get_kernel_name(nbfile)
        
        # run NB; trap when the kernel died (likely due to RAM exhaustion), or if kernel is not present
        try:
            nbpath, nb, errors = run_notebook(nbfile,kernel=kernel)
        except DeadKernelError:
            cprint('KERNEL DIED (LIKELY RAM EXHAUSTED)','red')
            test = 'FAIL'
        except NoSuchKernel:
            cprint('REQUIRED KERNEL NOT PRESENT','red')
            test = 'FAIL'
        except Exception as e:
            cprint('AN EXCEPTION OCCURRED: %s' % e,'red')
            test = 'FAIL'
        else:
            try:
                assert errors == []
            except:
                for e in errors:
                    traceback = e.pop('traceback')
                    for _ in traceback:
                        if plain is True:
                            _ = ansi_escape.sub('', _)
                        print(_)
                    print()

                cprint('TEST FAILURES ENCOUNTERED IN NOTEBOOK EXECUTION','red')
                test = 'FAIL'
            else:
                cprint('NOTEBOOK EXECUTED WITHOUT ERRORS','green')
                test = 'PASS'
        
        stop = time.time()
        duration = stop-start
        tests.append([nbfile,test,duration])  # append test results for this notebook
        print('RUNTIME: %.1f seconds' % duration)
        print('================================================')
        print('\n')
        
    # make a test result matrix
    df = pd.DataFrame(tests,columns=['notebook','test','duration'])
    subset = ['test']
    
    # styling
    bgcolor = lambda x: 'background-color: red' if x == 'FAIL' else 'background-color: green'
    fgcolor = lambda x: 'color: white'
    fontweight = lambda x: 'font-weight: bold'
    
    testresults = df.style.map(bgcolor,subset=subset).map(fgcolor,subset=subset).map(fontweight,subset=subset)
    
    return testresults


if __name__ == '__main__':

    start = time.time()
    
    cprint('Running testnotebooks.py',color='yellow',bar=5,newline=True)

    # List all paths with notebooks to test.
    # Paths can be single notebooks, directories, or a mix.
    # If single ipynb or single directory, paths can be a string.
    # If multiple elements, make it a sequence of strings.
    #
    # Examples:
    #
    #paths = '../01_GettingStartedWithDataLab/'  # test ony notebooks in 01_GettingStartedWithDataLab/ directory
    #paths = ('../01_GettingStartedWithDataLab/','../02_DataAccessOverview/')  # test notebooks in these two directories
    #paths = ('../01_GettingStartedWithDataLab/02_GettingStartedWithDataLab.ipynb','../02_DataAccessOverview/')  # test one NB in 01... dir, and all NBs in 02... dir
    paths = '../' # test all notebooks that are not excluded below

    # List pattern of notebook names, and/or paths, to exclude from testing.
    # ** means "any number of intermediate directories"
    # * means "any number of intermediate characters in this directory"
    # A few default notebooks are explicitly excluded here b/c they absolutely require interactive execution.
    # A few others are temporarily excluded until they will be fixed.
    exclude = ('**/*AuthClient.ipynb','**/Rowstore*.ipynb','**/e-Teen*/**/*ipynb',
               '**/DRAGONS_reduction_examples/GMOS_longslit_WhiteDwarf/GMOS_longslit_WhiteDwarf.ipynb',
               '**/DRAGONS_reduction_examples/GSAOI_Imaging_EllipticalGalaxy/GSAOI_Imaging_EllipticalGalaxy.ipynb')
    
    # log in to Data Lab once
    cprint('Login to Data Lab',color='yellow',bar=0,pad='')
    token = login(input('Username: '),getpass('Password: '))

    # get list of notebooks to run
    nbs = get_nbs(paths=paths,exclude=exclude)
    
    # run the tests (this can take a while)
    testresults = run(nbs,plain=True)

    stop = time.time()
    cprint('Total runtime: %g seconds' % (stop-start),color='yellow',bar=0,pad='',newline=True)
    
    cprint('TEST RESULT MATRIX',color='yellow',bar=0,pad='')
    print(testresults.data)
    print()

    cprint('All done.',color='yellow',bar=5,newline=True)
