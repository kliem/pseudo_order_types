# PseudoOrderTypes
Generate all pseudo order types for small point sets in the plane

We generate all pseudo order types according to

*O. Aichholzer, F. Aurenhammer, and H. Krasser. Enumerating Order Types for Small Point Sets with Applications. In Proc. 17th Ann. ACM Symp. Computational Geometry, pages 11-18, Medford, Massachusetts, USA, 2001.*

While the actual order types are listed here

http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/

the pseudo order types are not available anymore.

The package `pseudo_order_types` is `pip` installable and allows
iterating over pseudo order types with given number of points as well as
saving and loading them.

## Dependencies

- `setuptools`,
- `wheel`,
- `Cython`,
- `memory-allocator`.
