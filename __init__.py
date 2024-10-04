import os
import time
from typing import NoReturn, Sequence, Tuple, Optional
from dataclasses import dataclass, field
from contextlib import contextmanager
import traceback
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_fixed
import ase.db as db
import pandas as pd
from sqlite3 import OperationalError


@dataclass
class reaction:
    reactants: Sequence[Tuple[str, float]]
    products: Sequence[Tuple[str, float]]
    experimental_ref: float

    def toStr(self) -> str:
        return ' + '.join([f'{n:.2g}{smi if smi != "cid281" else "C|||O"}' for smi, n in self.reactants]) + ' ---> ' + ' + '.join([f'{n:.2g}{smi  if smi != "cid281" else "C|||O"}' for smi, n in self.products])


def build_pd(db_dir_list, select_key: Optional = None):
    if isinstance(db_dir_list, str): db_dir_list = [db_dir_list]
    for db_dir in db_dir_list:
        if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'):
            raise FileNotFoundError("Can't find database")
    db_list = [db.connect(work_db) for work_db in db_dir_list]
    pd_dat = pd.DataFrame([dict(atoms=row.toatoms(), **row.__dict__) for work_db in db_list for row in work_db.select(selection=select_key)])
    return pd_dat


def sanitize(unclean_str: str) -> str:
    for ch in ['!', '*', '?', '{', '[', '(', ')', ']', '}', "'", '"']: unclean_str = unclean_str.replace(ch, '')
    for ch in ['/', '\\', '|', ' ', ',', '.']: unclean_str = unclean_str.replace(ch, '_')
    for ch in ['=', '+', ':', ';']: unclean_str = unclean_str.replace(ch, '-')
    return unclean_str


@retry(retry=retry_if_exception_type(FileExistsError), stop=stop_after_attempt(5), wait=wait_fixed(2))
def folder_exist(folder_name: str, path: str = '.') -> None:
    if folder_name not in os.listdir(path): os.mkdir(ends_with(path, '/') + folder_name)


def ends_with(string: str, end_str: str) -> str:
    return string + end_str * (end_str != string[-len(end_str):0])


@retry(retry=retry_if_exception_type(OperationalError), stop=stop_after_attempt(5), wait=wait_fixed(10))
def update_db(db_dir: str, db_update_args: dict):
    with db.connect(db_dir) as db_obj:
        db_obj.update(**db_update_args)
