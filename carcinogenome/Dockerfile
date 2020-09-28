# Get shiny+tidyverse packages from rocker image
FROM rocker/shiny-verse:latest

# Set up the maintainer information
MAINTAINER Reina Chau (lilychau999@gmail.com)
    
# Set up a volume directory
VOLUME /srv/shiny-server/   

# Set up working directory to the app
WORKDIR /srv/shiny-server/

# Make the ShinyApp available at port 7856
EXPOSE 7856

# Copy configuration files to Docker image
COPY shiny-server-carcinogenome.sh /usr/bin/shiny-server.sh

# Allow permission
RUN ["chmod", "+rwx", "/srv/shiny-server/"]
RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]

CMD ["/usr/bin/shiny-server.sh"]

