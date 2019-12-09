SELECT * FROM EBI_PROJECT_DATA;
SELECT * FROM METADATA_EBI_PROJECT_DATA;
SELECT * FROM METADATA;
SELECT  * FROM METADATA_EBI_PROJECT_DATA;

-- kell egy subsampling-et csinálni. 

-- restarting subset info:
TRUNCATE TABLE SEQUENCING_SUBSET_MACHINE CASCADE;
TRUNCATE TABLE SEQUENCING_SUBSET_EBI_PROJECT_DATA CASCADE;
TRUNCATE TABLE sequencing_subset CASCADE;
TRUNCATE TABLE SUBSET_GROUP CASCADE;

DELETE FROM sequencing_subset PURGE;
DELETE FROM SUBSET_GROUP  PURGE;

COMMIT;




-- adding subsets 
INSERT INTO subset_group (base_name, description) VALUES ('TMP0', 'Test dataset');
commit;




SELECT * FROM subset_group;
SELECT * FROM tmp_sequencing_subsets;

--sequence_subset


drop table tmp_sequencing_subsets;
CREATE TABLE tmp_sequencing_subsets AS
SELECT  
21 subset_group_id,
'TMP0_' || T00.dataset_name || '_sscount' ||  floor(rownum/10) as subset_name, 
T00.EBI_PROJECT_DATA_ID
FROM (
SELECT 
t2.id EBI_PROJECT_DATA_ID,
t2.dataset_name
FROM METADATA t0
INNER JOIN METADATA_EBI_PROJECT_DATA t1
ON t1.METADATA_ID = t0.id
INNER JOIN EBI_PROJECT_DATA t2
ON t2.id = t1.ebi_project_data_id
WHERE t0.body_site = 'stool' 
and t0.dataset_id !=6
ORDER BY t2.dataset_id, t2.id
) T00
FETCH FIRST 100 ROWS ONLY;
;   

--truncate table SEQUENCING_SUBSET_EBI_PROJECT_DATA;
--drop table SEQUENCING_SUBSET_EBI_PROJECT_DATA;
--truncate table SEQUENCING_SUBSET;
--delete from SEQUENCING_SUBSET;

insert into SEQUENCING_SUBSET
(subset_group_id, subset_name)
SELECT distinct subset_group_id, subset_name
from tmp_sequencing_subsets;

INSERT INTO SEQUENCING_SUBSET_EBI_PROJECT_DATA
(SEQUENCING_SUBSET_ID, SUBSET_NAME,EBI_PROJECT_DATA_ID)
SELECT 
t0.id,
t0.SUBSET_NAME,
t1.EBI_PROJECT_DATA_ID
FROM SEQUENCING_SUBSET t0
INNER JOIN tmp_sequencing_subsets t1
ON t1.subset_group_id = t0.subset_group_id 
and t0.subset_name = t1.subset_name
;
DROP TABLE tmp_sequencing_subsets;

commit;

-- Adding stuff


