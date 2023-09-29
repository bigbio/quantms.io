The parquet format
------------------

Parquet is a column-based file format that can efficiently store and
process large amounts of data, and supports high compression and query
optimization.

Features
~~~~~~~~

-  **Column storage format**: Compared with JSON or CSV or any other row
   storage format, the main characteristic of Parquet is its design of
   column storage. This method can achieve more efficient compression
   and encoding, and improve the query performance of analytical
   statements. In general, tables in big data systems are stored as wide
   tables; that is, a table has many columns(dozens, even hundreds). The
   general query statements, especially analytical statements, only
   involve a small number of column queries. In this case, Parquet can
   read only the required columns, significantly reducing I/O and
   speeding up query execution.

   .. raw:: html

      <p align="center">

   .. raw:: html

      </p>

-  **Efficient compression and coding**: Parquet’s column storage format
   allows for better compression ratios because data in the same column
   tends to be more homogeneous. Parquet supports a variety of
   compression algorithms such as Snappy, Gzip, and LZO. In addition,
   Parquet uses advanced coding techniques such as RLE, bitpacking, and
   dictionary-encoding to further reduce storage requirements and
   improve query performance.

   .. raw:: html

      <p align="center">

   .. raw:: html

      </p>

-  **Schema evolutionary support**: Parquet aims to be able to handle
   changes in data schema over time, which is critical for big data
   systems. It supports schema evolution by allowing columns to be
   added, deleted, or modified without affecting existing data.

-  **Support for complex data types**: Parquet supports a rich set of
   data types, including nested and repeated structures, and data types
   such as arrays, maps, and structs. This feature allows you to build
   complex hierarchical data and store it efficiently in a compact
   binary format.

Theme
-----

For a large number of data result files, we aim to find a unified format
to store them, strictly control the data format and reduce storage
space.
