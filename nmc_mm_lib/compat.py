"""
compat.py: Python 2/3 compatibility module
@author: Kenneth Perrine
@contact: kperrine@utexas.edu
@organization: Network Modeling Center, Center for Transportation Research,
    Cockrell School of Engineering, The University of Texas at Austin 
@version: 1.0

@copyright: (C) 2016, The University of Texas at Austin
@license: GPL v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Borrowed from http://legacy.python.org/dev/peps/pep-0469/ 
try:
    dict.iteritems
except AttributeError:
    # Python 3
    def iterkeys(d):
        return iter(d.keys())
    def itervalues(d):
        return iter(d.values())
    def iteritems(d):
        return iter(d.items())
    def listkeys(d):
        return list(d.keys())
    def listvalues(d):
        return list(d.values())
    def listitems(d):
        return list(d.items())
else:
    # Python 2
    def iterkeys(d):
        return d.iterkeys()
    def itervalues(d):
        return d.itervalues()
    def iteritems(d):
        return d.iteritems()
    def listkeys(d):
        return d.keys()
    def listvalues(d):
        return d.values()
    def listitems(d):
        return d.items()
    