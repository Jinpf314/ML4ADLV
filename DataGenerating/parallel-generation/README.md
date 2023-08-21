# Parallel Generation
A set of auxilliary computer programs used to generate datasets 2 and 3 used in
section 5. These can use several CPU cores in parallel, and thus can generate
the datasets in few hours on a standard laptop PC.

# Requirements
- Sagemath
- MySQL

# Usage
Create a MySQL database, e.g. named adlv. In order to obtain dataset 2, run the following programs in parallel:

```
python3 dataset2.py A4 --random-seed 1 &
python3 dataset2.py A4 --random-seed 2 &
python3 dataset2.py A4 --random-seed 3 &
...
```

You can specify the database connection using the arguments --db-host, --db-user, --db-pass etc.

Once the programs have finished, you may run
```
python3 export.py A4 --db-tablename database2
```
to get a CSV file with the generated data and relevant features. The generation process for dataset 3 is analogous.

# PostgreSQL support

In order to connect to a PostgreSQL database, install the python module `psycopg2` and supply the following arguments to all program calls:

```
python3 [...] --db-pythonmodule psycopg2 --db-port 5432
```
