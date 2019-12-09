
-- creating schema for the basic metagenome data

--id	dataset_name	pubmed_id	balazs_comment	original_comments
truncate table metadata_EBI_project_data;
truncate table EBI_project_data;

truncate table sequencing_subset_machine;
truncate table sequencing_subset_EBI_project_data;
truncate table sequencing_subset;
truncate table subset_group;
truncate table LOG_MACHINE_PROCESSING;


drop SEQUENCE sequence_subset_seq;
drop SEQUENCE sequencing_subset_seq; 
drop SEQUENCE sequencing_subset_EBI_project_dat_seq; 
drop SEQUENCE sequencing_subset_machine_seq;
drop sequence subset_group_seq;
drop sequence log_machine_processing_seq;


DROP TABLE sequencing_subset_machine;
DROP TABLE sequencing_subset_EBI_project_data;
DROP TABLE sequencing_subset;
DROP TABLE subset_group;
drop table log_machine_processing;


DROP TABLE metagenome_dataset_meta_files;
DROP TABLE metadata;
drop table metagenome_dataset cascade constraint;
drop table EBI_project_data;
DROP TABLE  curated_metagenomes_table;
DROP TABLE metadata_EBI_project_data;
DROP TABLE log_machine_processing;


CREATE TABLE metagenome_dataset(
	id NUMBER,
	dataset_name  VARCHAR2(64),
	pubmed_id NUMBER,
	dataset_source  VARCHAR2(64),
	balazs_comment  VARCHAR2(4000),
	original_comments VARCHAR2(4000)
);

CREATE TABLE metagenome_dataset_meta_files(
id NUMBER,
dataset_id NUMBER,
dataset_name  VARCHAR2(64),
dataset_description_file	VARCHAR2(255),
file_category	VARCHAR2(64),
dataset_description_file_path  VARCHAR2(255)
);



CREATE TABLE EBI_project_data (
    id NUMBER,
    dataset_id  NUMBER,
    dataset_name  VARCHAR2(64),
    study_accession  VARCHAR2(64),
    secondary_study_accession  VARCHAR2(64),
    sample_accession  VARCHAR2(64),
    secondary_sample_accession  VARCHAR2(64),
    experiment_accession  VARCHAR2(64),
    run_accession  VARCHAR2(64),
    submission_accession  VARCHAR2(64),
    tax_id NUMBER,
    scientific_name  VARCHAR2(2000),
    instrument_platform  VARCHAR2(64),
    instrument_model  VARCHAR2(64),
    library_name  VARCHAR2(64),
    nominal_length NUMBER,
    library_layout  VARCHAR2(64),
    library_strategy  VARCHAR2(64),
    library_source  VARCHAR2(64),
    library_selection  VARCHAR2(64),
    read_count  NUMBER,
    base_count  NUMBER,
    center_name  VARCHAR2(64),
    first_public  VARCHAR2(64),
    last_updated  VARCHAR2(64),
    experiment_title  VARCHAR2(4000),
    study_title  VARCHAR2(4000),
    study_alias  VARCHAR2(4000),
    experiment_alias  VARCHAR2(4000),
    run_alias  VARCHAR2(256),
    fastq_bytes VARCHAR2(500),
    fastq_md5  VARCHAR2(1000),
    fastq_ftp  VARCHAR2(4000),
    fastq_aspera  VARCHAR2(4000),
    fastq_galaxy  VARCHAR2(4000),
    submitted_bytes  VARCHAR2(500),
    submitted_md5  VARCHAR2(256),
    submitted_ftp  VARCHAR2(4000),
    submitted_aspera  VARCHAR2(4000),
    submitted_galaxy  VARCHAR2(400),
    submitted_format  VARCHAR2(64),
    sra_bytes  VARCHAR2(64),
    sra_md5  VARCHAR2(64),
    sra_ftp  VARCHAR2(4000),
    sra_aspera  VARCHAR2(4000),
    sra_galaxy  VARCHAR2(4000),
    cram_index_ftp  VARCHAR2(64),
    cram_index_aspera  VARCHAR2(64),
    cram_index_galaxy  VARCHAR2(64),
    sample_alias  VARCHAR2(64),
    broker_name  VARCHAR2(64),
    sample_title  VARCHAR2(4000),
    nominal_sdev NUMBER,
    first_created  VARCHAR2(64)
);
--raw metadata table
CREATE TABLE metadata(
    id  number,
    dataset_id  number,
    dataset_name  VARCHAR2(64),
    study_condition  VARCHAR2(64),
    body_site  VARCHAR2(64),
    sampleID  VARCHAR2(64),
    subjectID  VARCHAR2(64),
    PMID  VARCHAR2(64),
    NCBI_accession  VARCHAR2(4000),
    non_westernized  VARCHAR2(64),
    country  VARCHAR2(64),
    disease  VARCHAR2(64),
    disease_subtype  VARCHAR2(64),
    gender  VARCHAR2(64),
    age  NUMBER,
    age_category  VARCHAR2(64),
    BMI  VARCHAR2(64),
    antibiotics_current_use  VARCHAR2(255),
    antibiotics_family  VARCHAR2(255),
    c_peptide  VARCHAR2(64),
    cholesterol  VARCHAR2(64),
    days_from_first_collection  VARCHAR2(64),
    lactating  VARCHAR2(64),
    pregnant  VARCHAR2(64),
    DNA_extraction_kit  VARCHAR2(64),
    sequencing_platform  VARCHAR2(64),
    median_read_length number,
    number_bases number,
    number_reads  number
);

CREATE TABLE metadata_EBI_project_data(
	metadata_id	number,        
    EBI_project_data_id NUMBER
);

-- Table for storing subsets for the batch processings
DROP TABLE sequencing_subset_machine;
DROP TABLE sequencing_subset_EBI_project_data;

CREATE TABLE subset_group
(
    id NUMBER,
    base_name VARCHAR(50),
    description VARCHAR(255),
    creation_date DATE
);

CREATE TABLE sequencing_subset(
	id NUMBER,
	subset_group_id NUMBER,
	subset_name VARCHAR(255),
	creation_date DATE
);


CREATE TABLE sequencing_subset_EBI_project_data(
	id NUMBER,
	sequencing_subset_id NUMBER,
	subset_name VARCHAR(255),
	EBI_project_data_id NUMBER,
	creation_date DATE
);

--mapping the subset to machines
CREATE TABLE sequencing_subset_machine(
 id NUMBER,
 sequencing_subset_id NUMBER,
 machine_name VARCHAR(255), --possibly the hostnames
 status VARCHAR(255),
 ASSIGNING_TIME DATE
);

drop table log_machine_processing;
CREATE TABLE log_machine_processing(
id NUMBER,
machine_name VARCHAR(255),
SEQUENCING_SUBSET_id NUMBER, -- ('This should be a valid subset id')
log_entry VARCHAR2(4000), --it could be very large I assume
log_date DATE
)
;

-- Adding constraints

ALTER TABLE metagenome_dataset
ADD CONSTRAINT metagenome_dataset_pk PRIMARY KEY (ID); 

ALTER TABLE metagenome_dataset
ADD CONSTRAINT metagenome_dataset_uq UNIQUE (dataset_name); 

ALTER TABLE metagenome_dataset
ADD CONSTRAINT metagenome_dataset_uq1 UNIQUE (id, dataset_name); 

ALTER TABLE metadata
ADD CONSTRAINT metadata_pk PRIMARY KEY (ID); 

ALTER TABLE metadata_EBI_project_data
ADD CONSTRAINT metadata_EBI_project_data_pk PRIMARY KEY (metadata_id, EBI_project_data_id); 

ALTER TABLE subset_group
ADD CONSTRAINT subset_group_pk PRIMARY KEY (id); 

ALTER TABLE sequencing_subset
ADD CONSTRAINT sequencing_subset_pk PRIMARY KEY (id); 

ALTER TABLE sequencing_subset_EBI_project_data
ADD CONSTRAINT sequencing_subset_EBI_project_data_pk PRIMARY KEY (id); 

ALTER TABLE sequencing_subset_machine
ADD CONSTRAINT sequencing_subset_machine_pk PRIMARY KEY (id); 

ALTER TABLE log_machine_processing
ADD CONSTRAINT log_machine_processing_pk PRIMARY KEY (id); 

ALTER TABLE metagenome_dataset_meta_files
ADD CONSTRAINT metagenome_dataset_meta_files_pk1 PRIMARY KEY (ID); 

-- EBI_project_data
--ALTER TABLE EBI_project_data
--ADD CONSTRAINT EBI_project_data_uq1 UNIQUE (run_accession); 

ALTER TABLE EBI_project_data
ADD CONSTRAINT EBI_project_data_pk PRIMARY KEY (id); 

ALTER TABLE EBI_project_data
ADD CONSTRAINT EBI_project_data_uq UNIQUE (run_accession); 

ALTER TABLE EBI_project_data
ADD CONSTRAINT EBI_project_data_uq2 UNIQUE (id,run_accession); 

ALTER TABLE sequencing_subset
ADD CONSTRAINT sequencing_subset_uq UNIQUE (subset_name); 

alter table EBI_project_data
add constraint EBI_project_data_run_acc_not_null check ( run_accession is not null );

alter table sequence_subset
add constraint sequence_subset_not_null check ( base_name is not null );

ALTER TABLE SEQUENCING_SUBSET_machine
ADD CONSTRAINT SEQUENCING_SUBSET_machines_uk unique (sequencing_subset_id);



ALTER TABLE metagenome_dataset_meta_files
ADD CONSTRAINT metagenome_dataset_meta_files_fk
  FOREIGN KEY (dataset_id)
  REFERENCES metagenome_dataset(id);


ALTER TABLE EBI_project_data
ADD CONSTRAINT EBI_project_data_fk1
 FOREIGN KEY (dataset_id)
REFERENCES metagenome_dataset(id);

ALTER TABLE metadata
ADD CONSTRAINT metadata_fk1
FOREIGN KEY (dataset_id)
REFERENCES metagenome_dataset(id);

ALTER TABLE metadata_EBI_project_data
ADD CONSTRAINT metadata_EBI_project_data_fk1
 FOREIGN KEY (EBI_project_data_id)
REFERENCES EBI_project_data(id);

ALTER TABLE metadata_EBI_project_data
ADD CONSTRAINT metadata_EBI_project_data_fk2
 FOREIGN KEY (metadata_id)
REFERENCES metadata(id);


ALTER TABLE sequencing_subset
ADD CONSTRAINT sequencing_subset_fk1
FOREIGN KEY (subset_group_id)
REFERENCES subset_group(id);

ALTER TABLE sequencing_subset_ebi_project_data
ADD CONSTRAINT sequencing_subset_prj_data_fk1
FOREIGN KEY (sequencing_subset_id)
REFERENCES sequencing_subset(id);

ALTER TABLE sequencing_subset_ebi_project_data
ADD CONSTRAINT sequencing_subset_prj_data_fk2
FOREIGN KEY (ebi_project_data_id)
REFERENCES ebi_project_data(id);

ALTER TABLE sequencing_subset_machine
ADD CONSTRAINT sequencing_subset_machine_fk1
FOREIGN KEY (sequencing_subset_id)
REFERENCES sequencing_subset(id);

ALTER TABLE log_machine_processing
ADD CONSTRAINT log_machine_processing_fk1
FOREIGN KEY (sequencing_subset_id)
REFERENCES sequencing_subset(id);

-- adding SEQUENCES:

CREATE SEQUENCE subset_group_seq START WITH 1; 
CREATE SEQUENCE sequencing_subset_seq START WITH 1; 
CREATE SEQUENCE sequencing_subset_EBI_project_dat_seq START WITH 1; 
CREATE SEQUENCE sequencing_subset_machine_seq START WITH 1; 
CREATE SEQUENCE log_machine_processing_seq START WITH 1; 






-- CREATING TRIGGERS --
/
CREATE OR REPLACE TRIGGER subset_group_TRIGGER
BEFORE INSERT ON subset_group
FOR EACH ROW
BEGIN
	IF :new.ID IS NULL THEN
		SELECT subset_group_seq.nextval INTO :new.id FROM DUAL;
	END IF;
END;
/
CREATE OR REPLACE TRIGGER subset_group_creation_date BEFORE INSERT ON subset_group 
FOR EACH ROW
BEGIN
 :NEW.creation_date := systimestamp;
END;
/
-- SEQUENCING SUBSETS --
/
CREATE OR REPLACE TRIGGER sequencing_subset_TRIGGER
BEFORE INSERT ON sequencing_subset
FOR EACH ROW
BEGIN
	IF :new.ID IS NULL THEN
		SELECT sequencing_subset_seq.nextval INTO :new.id FROM DUAL;
	END IF;
END;
/
CREATE OR REPLACE TRIGGER sequencing_subset_creation_date BEFORE INSERT ON sequencing_subset 
FOR EACH ROW
BEGIN
 :NEW.creation_date := systimestamp;
END;
/
CREATE OR REPLACE TRIGGER sequencing_subset_ebi_project_data_ID_TR
BEFORE INSERT ON SEQUENCING_SUBSET_EBI_PROJECT_DATA
FOR EACH ROW
BEGIN
	IF :new.ID IS NULL THEN
		SELECT sequencing_subset_EBI_project_dat_seq.nextval INTO :new.id FROM DUAL;
	END IF;
END;
/
CREATE OR REPLACE TRIGGER sequencing_subset_ebi_project_data_cr_date BEFORE INSERT ON sequencing_subset_ebi_project_data 
FOR EACH ROW
BEGIN
 :NEW.creation_date := systimestamp;
END;
/
/
CREATE OR REPLACE TRIGGER sequencing_subset_machine_ID_TR
BEFORE INSERT ON sequencing_subset_machine
FOR EACH ROW
BEGIN
	IF :new.ID IS NULL THEN
		SELECT sequencing_subset_machine_seq.nextval INTO :new.id FROM DUAL;
	END IF;
END;
/
CREATE OR REPLACE TRIGGER sequencing_subset_machine_cr_date BEFORE INSERT ON sequencing_subset_machine 
FOR EACH ROW
BEGIN
 :NEW.assigning_time := systimestamp;
END;
/
-- log_machine_processing
/
CREATE OR REPLACE TRIGGER log_machine_processing_ID_TR
BEFORE INSERT ON log_machine_processing
FOR EACH ROW
BEGIN
	IF :new.ID IS NULL THEN
		SELECT log_machine_processing_seq.nextval INTO :new.id FROM DUAL;
	END IF;
END;
/
CREATE OR REPLACE TRIGGER log_machine_processing_cr_date BEFORE INSERT ON log_machine_processing 
FOR EACH ROW
BEGIN
 :NEW.log_date := systimestamp;
END;
/



-- ADDING COMMENTS to tables

SELECT metadata_id, run_accession, count(*)
FROM metadata_EBI_project_data
GROUP BY metadata_id, run_accession
HAVING count(*)>1;


CREATE INDEX EBI_PROJECT_DATA_idx0 ON  EBI_PROJECT_DATA(SECONDARY_SAMPLE_ACCESSION);
CREATE BITMAP INDEX EBI_project_data_bmidx0 ON EBI_PROJECT_DATA (dataset_id);
CREATE BITMAP INDEX EBI_project_data_bmidx1 ON EBI_PROJECT_DATA (dataset_name);


CREATE INDEX METADATA_EBI_PROJECT_DATA_idx0 ON  METADATA_EBI_PROJECT_DATA(metadata_id);
CREATE INDEX METADATA_EBI_PROJECT_DATA_idx1 ON  METADATA_EBI_PROJECT_DATA(EBI_project_data_id);
CREATE INDEX METADATA_EBI_PROJECT_DATA_idx2 ON  METADATA_EBI_PROJECT_DATA(metadata_id,EBI_project_data_id);
CREATE INDEX METADATA_EBI_PROJECT_DATA_idx3 ON  METADATA_EBI_PROJECT_DATA(EBI_project_data_id,metadata_id);

CREATE BITMAP INDEX metadata_bmidx0 ON metadata (dataset_id);
CREATE BITMAP INDEX metadata_bmidx1 ON metadata (dataset_name);
CREATE BITMAP INDEX metadata_bmidx2 ON metadata (body_site);
CREATE BITMAP INDEX metadata_bmidx3 ON metadata (gender);

CREATE BITMAP INDEX metadata_dx0 ON metadata (body_site, id);


--CREATE BITMAP INDEX Sequencing_data_subsets_bmidx0 ON Sequencing_data_subsets (subset_base_name);
--CREATE BITMAP INDEX Sequencing_data_subsets_bmidx1 ON Sequencing_data_subsets (subset_base_name,sequence_subset_id);
--CREATE INDEX Sequencing_data_subsets_idx1 ON Sequencing_data_subsets (subset_base_name, id);


---- Procedures and functions
/
CREATE OR REPLACE PROCEDURE get_next_subset(query_subset_group_id IN NUMBER,subset_id OUT number) IS
CURSOR c IS 
   SELECT id
   FROM sequencing_subset ss  
   WHERE subset_group_id = query_subset_group_id 
    and   NOT EXISTS 
                  (SELECT id  
                   FROM SEQUENCING_SUBSET_MACHINE ssm 
                   WHERE ssm.sequencing_subset_id = ss.id )  
   for update skip locked;
free_id number;
BEGIN
OPEN c;
FOR i IN 1..1 LOOP
   FETCH c INTO free_id;
   EXIT WHEN c%NOTFOUND;
   select id into subset_id from sequencing_subset where id = free_id for update;
   DBMS_LOCK.sleep(0);
   
END LOOP;
CLOSE c;
END;
/
-- https://www.codeproject.com/Tips/778259/Using-an-Oracle-Database-Table-as-a-Multithreaded
;
/
exec dbms_lock.sleep(3); 
/
set serveroutput on;


/
CREATE OR REPLACE PROCEDURE Assign_datasubset_to_machine
  (machine_name IN VARCHAR,
  query_subset_group_id IN NUMBER,
  subset_id OUT number,
  new_association_id OUT number)
IS
  --new_association_id NUMBER := sequencing_subset_machine_seq.nextval;
BEGIN
  get_next_subset(query_subset_group_id,subset_id);
  IF subset_id IS NOT NULL then
       new_association_id :=sequencing_subset_machine_seq.nextval;
       INSERT INTO SEQUENCING_SUBSET_MACHINE (id, SEQUENCING_SUBSET_ID, machine_name, status)
       VALUES (new_association_id,subset_id,machine_name, 'STARTED');
       commit;
    END IF;
  --DBMS_OUTPUT.PUT_LINE(subset_id);
END;
/


DECLARE
  subset_id number;
  machine_name VARCHAR(255) := 'BELA';  
BEGIN
  get_next_subset(subset_id);
  INSERT INTO SEQUENCING_SUBSET_MACHINE (SEQUENCING_SUBSET_ID, machine_name)
  VALUES (subset_id,machine_name);
  commit;
  DBMS_OUTPUT.PUT_LINE(subset_id);
END; 
/


/
CREATE OR REPLACE FUNCTION Number_of_unprocessed_sequence_subsets 
( subset_group_id_to_check IN    NUMBER )
RETURN number IS 
   unprocessed_samples number := 0; 
BEGIN 
   SELECT
   count(distinct t0.id)- count(distinct t1.SEQUENCING_SUBSET_id) max_machines INTO unprocessed_samples
   FROM SEQUENCING_SUBSET t0
   LEFT OUTER JOIN SEQUENCING_SUBSET_MACHINE t1
   ON t1.SEQUENCING_SUBSET_id = t0.id
   WHERE subset_group_id = subset_group_id_to_check;
    
   RETURN unprocessed_samples; 
END; 
/ 


--- Building the query for the download
;
SELECT 
t1.subset_group_id,
t0.base_name,
t1.subset_name,
t1.id subset_id,
t3.run_accession sample_id,
t3.id EBI_PROJECT_DATA_ID ,
t3.dataset_id,
t3.dataset_name,
t3.library_strategy,
t3.fastq_aspera,
t3.fastq_md5
FROM subset_group t0
INNER JOIN SEQUENCING_SUBSET t1
ON t1.subset_group_id = t0.id
INNER JOIN SEQUENCING_SUBSET_EBI_PROJECT_DATA t2
ON t2.SEQUENCING_SUBSET_id = t1.id
INNER JOIN EBI_PROJECT_DATA t3
ON t3.id = t2.ebi_project_data_id
WHERE 
t1.id = 5045
and t0.id = 21
;

