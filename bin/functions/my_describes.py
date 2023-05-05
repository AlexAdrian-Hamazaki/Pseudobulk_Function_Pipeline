import numpy as np 


def describe(array):

    print(f'MEAN: {np.mean(array):.5f}')
    print(f'MEDIAN: {np.median(array):.5f}')
    print(f'MAX: {np.max(array):.5f}')
    print(f'MIN: {np.min(array):.5f}')

def describe_no_zeros(array):
    array = array[array!=0]
    print(f'MEAN: {np.mean(array):.5f}')
    print(f'MEDIAN: {np.median(array):.5f}')
    print(f'MAX: {np.max(array):.5f}')
    print(f'MIN: {np.min(array):.5f}')

def describe_no_nan(array):

    print(f'MEAN: {np.nanmean(array):.5f}')
    print(f'MEDIAN: {np.nanmedian(array):.5f}')
    print(f'MAX: {np.max(array):.5f}')
    print(f'MIN: {np.min(array):.5f}')