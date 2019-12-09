#! /bin/bash

echo "First arg: $1"
echo "Second arg: $2"

username=$1
tablespace=$2
defaultpasswd="BioProject"

echo "CREATE USER $username  PROFILE \"DEFAULT\" IDENTIFIED BY $defaultpasswd  DEFAULT TABLESPACE $tablespace TEMPORARY TABLESPACE \"TEMP\" ACCOUNT UNLOCK;"
echo "GRANT ALTER ANY INDEX TO $username;"  ;
echo "GRANT ALTER ANY INDEXTYPE TO $username;"  ;
echo "GRANT ALTER ANY PROCEDURE TO $username;"  ;
echo "GRANT ALTER ANY ROLE TO $username;"  ;
echo "GRANT ALTER ANY SEQUENCE TO $username;"  ;
echo "GRANT ALTER ANY TABLE TO $username;"  ;
echo "GRANT ANALYZE ANY TO $username;"  ;
echo "GRANT CREATE ANY INDEX TO $username;"  ;
echo "GRANT CREATE ANY INDEXTYPE TO $username;"  ;
echo "GRANT CREATE ANY PROCEDURE TO $username;"  ;
echo "GRANT CREATE ANY SEQUENCE TO $username;"  ;
echo "GRANT CREATE ANY SYNONYM TO $username;"  ;
echo "GRANT CREATE ANY TABLE TO $username;"  ;
echo "GRANT CREATE ANY TRIGGER TO $username;"  ;
echo "GRANT CREATE ANY TYPE TO $username;"  ;
echo "GRANT CREATE ANY VIEW TO $username;"  ;
echo "GRANT CREATE PROCEDURE TO $username;"  ;
echo "GRANT CREATE PUBLIC SYNONYM TO $username;"  ;
echo "GRANT CREATE ROLE TO $username;"  ;
echo "GRANT CREATE TABLE TO $username;"  ;
echo "GRANT CREATE VIEW TO $username;"  ;
echo "GRANT \"CONNECT\" TO $username;"  ;
echo "GRANT CREATE TYPE TO $username;"  ;
echo "GRANT CREATE ANY INDEX TO $username;"  ;
#echo "ALTER USER $username  QUOTA 10M ON USERS;";

