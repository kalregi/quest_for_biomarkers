a-- Adding a limited subset for testing purposes

--NAME (tiny test 1): TT1

COMMIT;

INSERT INTO subset_group (
    base_name,
    description
) VALUES (
    'TT1',
    'Small test set (1 set consists of 2 element)'
);

select * from subset_group;


COMMIT;

SELECT
    *
FROM
    subset_group;
    

DROP TABLE tmp_sequencing_subsets;

CREATE TABLE tmp_sequencing_subsets
    AS
        SELECT
            2 subset_group_id,
            'TT1_'
            || t00.dataset_name
            || '_sscount'
            || floor(ROWNUM / 2) AS subset_name,
            t00.ebi_project_data_id
        FROM
            (
                SELECT
                    t2.id ebi_project_data_id,
                    t2.dataset_name
                FROM
                    metadata t0
                    INNER JOIN metadata_ebi_project_data t1 ON t1.metadata_id = t0.id
                    INNER JOIN ebi_project_data t2 ON t2.id = t1.ebi_project_data_id
                WHERE
                    t0.body_site = 'stool'
                    AND   t0.dataset_id != 6
                ORDER BY
                    t2.dataset_id,
                    t2.id
            ) t00
        FETCH FIRST 6 ROWS
only;
;
--association table
SELECT * FROM tmp_sequencing_subsets;



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
-- 
SELECT * FROM SUBSET_GROUP;
SELECT * FROM SEQUENCING_SUBSET_machine;
SELECT * FROM SEQUENCING_SUBSET where subset_group_id = 2;

SELECT * FROM LOG_machine_processing ORDER BY ID DESC;
WHERE log_entry LIKE '%ERR209888%';





