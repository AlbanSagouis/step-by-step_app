#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)
# library(mobsim)


# Define UI for slider demo application
  
navbarPage("Visualization of biodiversity pattern", selected="Step-by-step",
	tabPanel("Step-by-step",
		sidebarLayout(
			sidebarPanel(
				fluidRow(
					column(width=6,
						# number of species
						numericInput("sbsS", "Species Richness", min=5, max=500, value=50, step=5)
					),
					column(width=6,
						# number of individuals
						numericInput("sbsN", "Number of individuals", min=10, max=5000, value=1000, step=10)
					)
				),
				selectizeInput("sbssad_type", "SAD Type", choices=c("lognormal"="lnorm","geometric"="geom","Fisher's log-series"="ls")),
				uiOutput("sbsCVslider"),

				selectizeInput(inputId="sbsspatdist", "Cluster parameter", choices = c("Number of mother points"="n.mother", "Number of clusters"="n.cluster")),
				helpText("Number of mother points per species OR number of individuals per cluster."),
				textInput(inputId="sbsspatcoef",label="Integer values separated by commas", value="0"),
				textInput(inputId="sbsspatagg", label="Spatial Aggregation (mean distance from mother points)", value = 0.1),

				
				# Restart action button
				actionButton(inputId="sbsRestart",label="Restart Simulation"),
				
				# sampling parameters
				fluidRow(
					column(width=6,
						numericInput("sbsnumber_of_quadrats", label="Number of quadrats", value=20, min=1, max=1000, step=1)
					),
					column(width=6,
						numericInput("sbsarea_of_quadrats", label="Area of quadrats", value=0.005, min=0.00001, max=1, step=0.005)
					)
				),
				actionButton("sbsnew_sampling_button", label="Restart sampling"),
				# Next step
				actionButton("sbskeep_step", label="Next step")
			),
			mainPanel(
			# plot of the community
			# distance decay
				uiOutput("sbsfirst_step")
			)
		)
	)
)