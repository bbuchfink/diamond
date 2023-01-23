"""
A python wrapper module of diamond
"""

from typing import Any
from .libdiamond import main, version
import os
import functools
import inspect


def _require_db(func):
    """
    An internal decorator for diamond to check database.
    """
    @functools.wraps(func)
    def wrapper_func(*args, **kwargs):
        if len(args) == 0 or not isinstance(args[0], Diamond):
            raise ValueError(
                "This decorator is only used for methods of `Diamond`")
        self: Diamond = args[0]
        if not self._check_db():
            raise RuntimeError("Database not existed, please check it first")
        return func(*args, **kwargs)
    return wrapper_func


def not_null(func):
    """
    A decorator for checking params.
    When parameter type is annotated, and default
    value is not `None` if given, this parameter is 
    required as annotation type.

    When type is not annotated or annotated as `Any`, or 
    default value is set to `None, we will allow it to
    be `None`.
    """
    @functools.wraps(func)
    def wrapper_func(*args, **kwargs):
        params = inspect.signature(func).parameters
        name_list = func.__code__.co_varnames

        def check_annotation(param_name, param_val):
            annotation_type: type = params[param_name].annotation
            default_value = params[param_name].default
            if annotation_type == inspect._empty or annotation_type == Any:
                return
            if default_value is None and param_val is None:
                return
            if str(annotation_type).startswith("typing."):
                if param_val is None:
                    raise ValueError(
                        f"Parameter '{param_name}' requires not be `NoneType`")

            elif not isinstance(param_val, annotation_type):
                raise ValueError(
                    f"Parameter '{param_name}' requires be an instance of {annotation_type}")

        # check position params
        for i, param_val in enumerate(args):
            param_name = name_list[i]
            check_annotation(param_name, param_val)

        # check keyword params
        for param_name, param_val in kwargs.items():
            check_annotation(param_name, param_val)
        return func(*args, **kwargs)
    return wrapper_func


class Diamond(object):
    """
    Diamond python wrapper object
    """

    @not_null
    def __init__(self, database: str, n_threads: int = 1,
                 quiet: bool = False, log: bool = False) -> None:

        self.database = database
        self.n_threads = n_threads
        self.quiet = quiet
        self.log = log

    @_require_db
    @not_null
    def blastp(self, query: str, out: str) -> int:
        """
        Sequence alignment for proteins by diamond.

        Parameters:
        ------------
        db: str
            the database file
        query: str
            the input fasta query file

        """
        args = self._add_general_options(
            "blastp", "--query", query, "--out", out)
        return main(*args)

    @_require_db
    @not_null
    def dbinfo(self):
        """
        Print the information of current database.
        """
        args = self._add_general_options("dbinfo")
        return main(*args)

    @not_null
    def makedb(self, input_seqs_file: str,
               taxonmap: str = None,
               taxonnodes: str = None,
               taxonnames: str = None):
        """
        Make database based on given FASTA files.

        Parameters:
        -------------
        input_seqs_file: str
            Path to the input protein reference database
            file in FASTA format (may be gzip compressed).
        
        taxonmap: str
            Path to mapping file that maps NCBI protein
            accession numbers to taxon ids (gzip compressed).
        
        taxonnodes: str
            Path to the nodes.dmp file from the NCBI taxonomy.
        
        taxonnames: str
            Path to the names.dmp file from the NCBI taxonomy. 
        """
        args = self._add_general_options(
            "makedb", "--in", input_seqs_file,
            taxonmap=taxonmap,
            taxonnodes=taxonnodes,
            taxonnames=taxonnames)
        return main(*args)

    def _add_general_options(self, *args, **kwargs) -> list:
        new_args = list(args)
        kwargs["db"] = self.database
        kwargs["threads"] = str(self.n_threads)
        for k, v in kwargs.items():
            if v:
                new_args.append(f"--{k}")
                new_args.append(v)

        if self.quiet:
            new_args.append("--quiet")
        if self.log:
            new_args.append("--log")
        return new_args

    def _check_db(self) -> bool:
        return os.path.exists(self.database)

    @property
    def version(self) -> str:
        """Version of diamond"""
        return version()
