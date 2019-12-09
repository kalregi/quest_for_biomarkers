---- SCHEMA FOR MANAGING google-cloud VM-s
-- juhjkhkjh
-- cretate vm data
-- hghjgj
-- jkhkjh
-- ljhkjh

-- google cloud auto shutdow feature.

truncate table machine_status;
truncate table machine;
truncate table machine_group;

drop table machine_status;
drop table machine;
DROP TABLE machine_group;

drop sequence machine_group_seq;
drop sequence machine_seq;
drop sequence machine_status_seq;





CREATE TABLE machine_group(
id NUMBER,
short_name VARCHAR(8), --maximum 8char prefix for the VM-s belonging to that group
DESCRIPTION varchar(255), -- DESCRIPTION ABOUT the vm group
creation_date DATE
);

CREATE TABLE machine(
 id NUMBER,
 machine_group_id NUMBER,
 host_name VARCHAR(255) NOT NULL, --short hostname
 machine_prefix VARCHAR(255), --short prefix for the group name
 machine_name_analysis VARCHAR(255), --machine name used in the analysis (it could be a random number, normally it should contain the id), should be unique
 cloud_zone VARCHAR(255),
 machine_type VARCHAR(255), --type of the machine
 preemptible VARCHAR(255), --preemptibility of the machine
 internal_ip VARCHAR(16), --internal IP address of the machine
 external_ip VARCHAR(16), --external IP
 creation_date  DATE, --first date of the registration
 act_status VARCHAR(255) --actual state
  );

CREATE TABLE machine_status(
 id NUMBER,
 machine_id NUMBER,
 machine_status VARCHAR(255), --logging the status change
 status_registration DATE
);

--adding constraints 
ALTER TABLE machine_group
ADD CONSTRAINT machine_group_pk PRIMARY KEY (ID); 
ALTER TABLE machine
ADD CONSTRAINT machine_pk PRIMARY KEY (ID); 
ALTER TABLE machine_status
ADD CONSTRAINT machine_status_pk PRIMARY KEY (ID); 
-- 
ALTER TABLE machine
ADD CONSTRAINT machine_uq UNIQUE (machine_name_analysis); 
-- FOREIGN keys

ALTER TABLE machine
ADD CONSTRAINT machine_fk
FOREIGN KEY (machine_group_id)
REFERENCES machine_group(id);

ALTER TABLE machine_status
ADD CONSTRAINT machine_status_fk
FOREIGN KEY (machine_id)
REFERENCES machine(id);

-- sequences --

CREATE SEQUENCE machine_group_seq START WITH 1; 
CREATE SEQUENCE machine_seq START WITH 1; 
CREATE SEQUENCE machine_status_seq START WITH 1; 

-- ADDING TRIGGERS  --




--
SELECT
id,
machine_name,
log_date,
TO_CHAR(log_date, 'YYYY-MON-DD HH24:MI:SS')

FROM LOG_MACHINE_PROCESSING
;

