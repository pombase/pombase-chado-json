--
-- PostgreSQL database dump
--

-- Dumped from database version 12.2 (Ubuntu 12.2-2.pgdg19.10+1)
-- Dumped by pg_dump version 12.2 (Ubuntu 12.2-2.pgdg19.10+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: pgcrypto; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pgcrypto WITH SCHEMA public;


--
-- Name: EXTENSION pgcrypto; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION pgcrypto IS 'cryptographic functions';


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: query; Type: TABLE; Schema: public; Owner: kmr44
--

CREATE TABLE public.query (
    id uuid NOT NULL,
    query_json jsonb NOT NULL,
    creation_timestamp timestamp without time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.query OWNER TO kmr44;

--
-- Name: query query_id_key; Type: CONSTRAINT; Schema: public; Owner: kmr44
--

ALTER TABLE ONLY public.query
    ADD CONSTRAINT query_id_key UNIQUE (id);


--
-- Name: query_json_sha256_idx; Type: INDEX; Schema: public; Owner: kmr44
--

CREATE UNIQUE INDEX query_json_sha256_idx ON public.query USING btree (public.digest((query_json)::text, 'sha256'::text));


--
-- PostgreSQL database dump complete
--

