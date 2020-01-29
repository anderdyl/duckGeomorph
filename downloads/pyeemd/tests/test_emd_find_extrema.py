
# Copyright 2013 Perttu Luukko

# This file is part of libeemd.
# 
# libeemd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# libeemd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with libeemd.  If not, see <http://www.gnu.org/licenses/>.

from pyeemd import emd_find_extrema
from nose.tools import assert_equal, assert_true, assert_false, assert_greater, assert_false, raises
from numpy import zeros, all, repeat
from numpy.random import normal, randint

@raises(ValueError)
def test_bogus():
    x = emd_find_extrema("This spoon is too big")

@raises(ValueError)
def test_wrong_dims():
    x = zeros((2,2))
    emd_find_extrema(x)

def test_empty():
    x = []
    maxx, maxy, minx, miny = emd_find_extrema(x)
    assert_equal(maxx.size, 0)
    assert_equal(maxy.size, 0)
    assert_equal(minx.size, 0)
    assert_equal(miny.size, 0)

def test_size_1():
    x = [2]
    maxx, maxy, minx, miny = emd_find_extrema(x)
    assert all(maxx == [0])
    assert all(maxy == [2])
    assert all(minx == [0])
    assert all(miny == [2])

def test_size_2():
    x = [2, 5]
    maxx, maxy, minx, miny = emd_find_extrema(x)
    assert all(maxx == [0, 1])
    assert all(maxy == x)
    assert all(minx == [0, 1])
    assert all(miny == x)

def test_size_3():
    x = [2, -1, 5]
    maxx, maxy, minx, miny = emd_find_extrema(x)
    assert all(maxx == [0, 2])
    assert all(maxy == [2, 5])
    assert all(minx == [0, 1, 2])
    assert all(miny == x)

def test_zeros():
    zs = zeros(10)
    maxx, maxy, minx, miny = emd_find_extrema(zs)
    assert all(maxx == [0, 9])
    assert all(maxy == [0, 0])
    assert all(minx == [0, 9])
    assert all(miny == [0, 0])

def test_extrema():
    for i in range(16):
        x = normal(0, 1, 128)
        maxx, maxy, minx, miny = emd_find_extrema(x)
        yield check_extrema, x, maxx, maxy, minx, miny

def test_extrema_flats():
    for i in range(16):
        x = repeat(normal(0, 1, 64), randint(1, 3, 64))
        maxx, maxy, minx, miny = emd_find_extrema(x)
        yield check_extrema, x, maxx, maxy, minx, miny

def check_extrema(x, maxx, maxy, minx, miny):
    assert_equal(maxx[0], 0)
    assert_equal(maxx[-1], len(x)-1)
    assert_equal(minx[0], 0)
    assert_equal(minx[-1], len(x)-1)
    for n, i in enumerate(maxx):
        if (i == 0 or i == len(x)-1):
            continue
        if (i == int(i)):
            i = int(i)
            assert_equal(maxy[n], x[i])
            assert x[i-1] < x[i]
            assert x[i+1] < x[i]
        else:
            i = int(i)
            assert_equal(x[i], x[i])
            assert_equal(x[i+1], x[i])
    for n, i in enumerate(minx):
        if (i == 0 or i == len(x)-1):
            continue
        if (i == int(i)):
            i = int(i)
            assert_equal(miny[n], x[i])
            assert x[i-1] > x[i]
            assert x[i+1] > x[i]
        else:
            i = int(i)
            assert_equal(x[i], x[i])
            assert_equal(x[i+1], x[i])
