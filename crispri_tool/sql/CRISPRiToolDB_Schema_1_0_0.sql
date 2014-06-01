
CREATE DATABASE IF NOT EXISTS CRISPRiToolDB;
USE CRISPRiToolDB;

CREATE TABLE IF NOT EXISTS GeneDefinition
(
    GeneId          INT     NOT NULL PRIMARY KEY AUTO_INCREMENT,
    GeneSequence    TEXT    NOT NULL
);

CREATE TABLE IF NOT EXISTS StrandDefinition
(
    StrandId        INT     NOT NULL PRIMARY KEY,
    StrandName      TEXT    NOT NULL
);

CREATE TABLE IF NOT EXISTS GeneTarget
(
    TargetId        INT     NOT NULL PRIMARY KEY AUTO_INCREMENT,
    GeneId          INT     NOT NULL,
    StrandId        INT     NOT NULL,
    Position        INT     NOT NULL,
    Effectiveness   FLOAT   NULL,
    sgRNA           TEXT    NOT NULL,

    FOREIGN KEY( GeneId ) REFERENCES GeneDefinition( GeneId ),
    FOREIGN KEY( StrandId ) REFERENCES StrandDefinition( StrandId )
);

INSERT IGNORE INTO StrandDefinition VALUES
    ( 0, 'Coding' ),
    ( 1, 'Template' );
