# isomiRTar
## How to work with dockerized DB
First, cd to db/postgresql_docker. Then:
- docker-compose up -d (run docker image)
- docker-compose stop (stop container, but do not remove it)
- docker-compose down (stop container and remove it)
- docker ps (see running containers)
- docker exec -it containername /bin/bash (attach to container's shell)
- psql -h localhost -d dbname -U username (connect to the db)
