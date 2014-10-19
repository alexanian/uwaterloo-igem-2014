DROP DATABASE IF EXISTS possible_targets;

CREATE DATABASE possible_targets;
USE possible_targets;

CREATE TABLE GeneTable
(
    GeneId      INT             NOT NULL PRIMARY KEY AUTO_INCREMENT,
    GeneName    VARCHAR(512)    NOT NULL,
    CodingSeq   TEXT            NOT NULL,
    
    UNIQUE( GeneName )
);

CREATE TABLE StrandTable
(
    StrandId    INT             NOT NULL PRIMARY KEY AUTO_INCREMENT,
    Strand      VARCHAR(10)     NOT NULL
);

INSERT INTO StrandTable( Strand ) VALUES 
    ( 'Coding' ),
    ( 'Template' );

CREATE TABLE Targets
(
    GeneId      INT             NOT NULL,
    StrandId    INT             NOT NULL,
    Position    INT             NOT NULL,
    Result      TINYTEXT        NOT NULL,
    Repression  DOUBLE          NOT NULL DEFAULT 0,

    UNIQUE ( GeneId, StrandId, Position ),
    FOREIGN KEY ( GeneId ) REFERENCES GeneTable( GeneId ),
    FOREIGN KEY ( StrandId ) REFERENCES StrandTable( StrandId )
);
