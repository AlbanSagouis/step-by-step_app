
library(shiny)
# library(devtools)
# install_github('MoBiodiv/mobsim')    # downloads the latest version of the package
# library(mobsim, lib.loc="./Library")
source("code/Diversity_Area_Relationships.R", local = TRUE)
source("code/rThomas_r.r", local = TRUE)
source("code/Sample_quadrats.R", local = TRUE)
source("code/Sim_Community.R", local = TRUE)
# library(ggplot2)
# library(markdown)

# Define server logic for slider examples
shinyServer(function(input, output, session) {
  
	# STEP BY STEP
						
	## Parameters
	# update range for species richness, an observed species has minimum one individual
	observe({
		updateSliderInput(session, "sbsS", min=5, max=input$sbsN, value=50, step=5)
	})
	output$sbsCVslider <- renderUI({
		switch(input$sbssad_type,
			"lnorm"=sliderInput("sbscoef", label="CV(abundance), i.e. standard deviation of abundances divided by the mean abundance",value=1, min=0, max=5, step=0.1, ticks=F),
			"geom"=sliderInput("sbscoef", label="Probability of success in each trial. 0 < prob <= 1", value=0.5, min=0, max=1, step=0.1, ticks=F),
			"ls"=textInput("sbscoef", label="Fisher's alpha parameter", value=1)
		)
	})
	
	
	
	## community simulation
	sbssim.com <- reactive({
		input$sbsRestart
    
		isolate({
			# set.seed(229377)	# 229376
		
			spatagg_num <- as.numeric(unlist(strsplit(trimws(input$sbsspatagg), ",")))
			spatcoef_num <- as.numeric(unlist(strsplit(trimws(input$sbsspatcoef), ",")))
		 
			if(input$sbsspatdist=="n.mother") n.mother <- spatcoef_num else n.mother <- NA
			if(input$sbsspatdist=="n.cluster") n.cluster <- spatcoef_num else n.cluster <- NA
			
			simulation_parameters <- list(mother_points=n.mother,
																		cluster_points=n.cluster,
																		xmother=NA,
																		ymother=NA)
									

			switch(input$sbssad_type,
							"lnorm"=sim_thomas_community(s_pool = input$sbsS, n_sim = input$sbsN, 
								sigma=spatagg_num, mother_points=simulation_parameters$mother_points, cluster_points=simulation_parameters$cluster_points, xmother=simulation_parameters$xmother, ymother=simulation_parameters$ymother,
								sad_type = input$sbssad_type, sad_coef=list(cv_abund=input$sbscoef),
								fix_s_sim = T),
							"geom"=sim_thomas_community(s_pool = input$sbsS, n_sim = input$sbsN,
								sigma=spatagg_num, mother_points=simulation_parameters$mother_points, cluster_points=simulation_parameters$cluster_points, xmother=simulation_parameters$xmother, ymother=simulation_parameters$ymother,
								sad_type = input$sbssad_type, sad_coef=list(prob=input$sbscoef),
								fix_s_sim = T),
							"ls"=sim_thomas_community(s_pool = input$sbsS, n_sim = input$sbsN,
								sad_type = input$sbssad_type, sad_coef=list(N=input$sbsN,alpha=as.numeric(input$sbscoef)),
								sigma=spatagg_num, mother_points=simulation_parameters$mother_points, cluster_points=simulation_parameters$cluster_points, xmother=simulation_parameters$xmother, ymother=simulation_parameters$ymother,
								fix_s_sim = T)
						)

		# session$userData$sbsprevious.sim.com <- sbssim.com()
		
		})
	})

	## Community summary
	output$sbscommunity_summary_table <- renderTable({
		input$sbsRestart
		# data.frame(sbsspatagg=input$sbsspatagg,
			# sbscoef=input$sbscoef,
			# sbsspatcoef=input$sbsspatcoef,
			# sbsspatdist=input$sbsspatdist,
			# sbssad_type=input$sbssad_type,
			# sbsS=input$sbsS,
			# sbsN=input$sbsN)
		# data.frame(isnull=is.null(sbssim.com()),
			# class=class(sbssim.com()),
			# length=length(sbssim.com()))
		# sbssim.com()$census
		data.frame(Community = "",
						n_species = length(levels(sbssim.com()$census$species)),
						n_individuals = nrow(sbssim.com()$census))
	})
	
	## Sampling 
	sbssampling_quadrats <- reactive({
		input$sbsnew_sampling_button
		
		sample_quadrats(comm=sbssim.com(), n_quadrats=input$sbsnumber_of_quadrats, quadrat_area=input$sbsarea_of_quadrats, avoid_overlap=T, plot=F)
	})
	
	## Sampling summary
	### gamma scale
		output$sbsgamma_table <- renderTable({
		input$sbsRestart
		input$sbsnew_sampling_button
		
		isolate({
			abund <- apply(sbssampling_quadrats()$spec_dat, 2, sum)
			abund <- abund[abund > 0]
			relabund <- abund/sum(abund)
			shannon <- - sum(relabund * log(relabund))
			simpson <- 1- sum(relabund^2)
			data.frame(
							Gamma = "",
							n_species= sum(abund >0),
							shannon = round(shannon, 3),
							# ens_shannon = round(exp(shannon), 3),
							simpson = round(simpson, 3)
							# ens_simpson = round(1/(1 - simpson), 3)
			)
		})
	})
							

	### alpha scale
	output$sbsalpha_summary_table <- renderTable({
		input$sbsRestart
		input$sbsnew_sampling_button
		
		isolate({
			quadrats_coordinates <- sbssampling_quadrats()$xy_dat
			temp <- 	as.data.frame(t(round(sapply(1:nrow(quadrats_coordinates), function(i) {
				div_rect(x0=quadrats_coordinates$x[i], y0=quadrats_coordinates$y[i], xsize=sqrt(input$sbsarea_of_quadrats), ysize=sqrt(input$sbsarea_of_quadrats), comm=sbssim.com())
			}), 3)))[,c('n_species','n_endemics','shannon','simpson')]
			funs <- list(min=min, max=max, mean=mean, sd=sd)
			data.frame(Alpha=colnames(temp), round(sapply(funs, mapply, temp),3))
		})
	})

	
	## Plot the community and sampling squares
	
	output$sbssampling_plot <- renderPlot({
		# isolate({
			quadrats_coordinates <- sbssampling_quadrats()$xy_dat
			plot(sbssim.com(), main = "Community distribution")
			graphics::rect(quadrats_coordinates$x,
								quadrats_coordinates$y,
								quadrats_coordinates$x + sqrt(input$sbsarea_of_quadrats),
								quadrats_coordinates$y + sqrt(input$sbsarea_of_quadrats),
								lwd = 2, col = grDevices::adjustcolor("white", alpha.f = 0.6))
		# })
	})
	
	
	output$sbsdistance_decay_plot <- renderPlot({
		input$sbsRestart
		input$sbsnew_sampling_button
		
		isolate({
			plot(dist_decay_quadrats(sbssampling_quadrats(), method = "bray", binary = F))
		})
	})
	
	
	output$sbsfirst_step <- renderUI({
		fluidRow(align="center",
			column(width=4,
				tableOutput("sbscommunity_summary_table"),
				tableOutput("sbsgamma_table"),
				tableOutput("sbsalpha_summary_table")
			),
			column(width=4,
				plotOutput("sbssampling_plot")
			),
			column(width=4,
				plotOutput("sbsdistance_decay_plot")
			),
			hr()
		)
	})
	
	
	########## END OF FIRST STEP
	
	
	
	
	
	
	
	
	
	
	
	sbsvalues <- reactiveValues()
	
	observeEvent(input$sbskeep_step, {
		# ui ID
		divID <- gsub("\\.", "", format(Sys.time(), "%H%M%OS3"))
		btnID <- paste0(divID, "rmv")
		community_summary_id <- paste0(divID, "comm")
		gamma_table_id <- paste0(divID, "gamm")
		alpha_summary_table_id <- paste0(divID, "alph")
		sampling_plot_ID <- paste0(divID, "samp")
		distance_decay_plot_ID <- paste0(divID, "dist")
	
		# only create button if there is none
		if (is.null(sbsvalues[[divID]])) {
      
			insertUI(
				selector = "#sbsfirst_step",
				where="afterEnd",
				ui = tags$div(id = divID,
					fluidRow(
						column(width=4,
							tableOutput(community_summary_id),
							tableOutput(gamma_table_id),
							tableOutput(alpha_summary_table_id),
							actionButton(btnID, "Remove this step", class = "pull-right btn btn-danger")
						),
						column(width=4,
							plotOutput(sampling_plot_ID)
						),
						column(width=4,
							plotOutput(distance_decay_plot_ID)
						)
					)
				)
			)
      
      output[[community_summary_id]] <- renderTable({
			isolate({
				data.frame(Community = "",
							n_species = length(levels(sbssim.com()$census$species)),
							n_individuals = nrow(sbssim.com()$census))
			})
		})
		
      output[[gamma_table_id]] <- renderTable({
			isolate({
				abund <- apply(sbssampling_quadrats()$spec_dat, 2, sum)
				abund <- abund[abund > 0]
				relabund <- abund/sum(abund)
				shannon <- - sum(relabund * log(relabund))
				simpson <- 1- sum(relabund^2)
				data.frame(
								Gamma = "",
								n_species= sum(abund >0),
								shannon = round(shannon, 3),
								# ens_shannon = round(exp(shannon), 3),
								simpson = round(simpson, 3)
								# ens_simpson = round(1/(1 - simpson), 3)
				)
			})
		})      
		
		output[[alpha_summary_table_id]] <- renderTable({
			isolate({
				quadrats_coordinates <- sbssampling_quadrats()$xy_dat
				temp <- 	as.data.frame(t(round(sapply(1:nrow(quadrats_coordinates), function(i) {
					div_rect(x0=quadrats_coordinates$x[i], y0=quadrats_coordinates$y[i], xsize=sqrt(input$sbsarea_of_quadrats), ysize=sqrt(input$sbsarea_of_quadrats), comm=sbssim.com())
				}), 3)))[,c('n_species','n_endemics','shannon','simpson')]
				funs <- list(min=min, max=max, mean=mean, sd=sd)
				data.frame(Alpha=colnames(temp), round(sapply(funs, mapply, temp),3))
			})
		})
		
		output[[sampling_plot_ID]] <- renderPlot({
			isolate({
				quadrats_coordinates <- sbssampling_quadrats()$xy_dat
				plot(sbssim.com(), main = "Community distribution")
				graphics::rect(quadrats_coordinates$x,
									quadrats_coordinates$y,
									quadrats_coordinates$x + sqrt(input$sbsarea_of_quadrats),
									quadrats_coordinates$y + sqrt(input$sbsarea_of_quadrats),
									lwd = 2, col = grDevices::adjustcolor("white", alpha.f = 0.6))
			})
		})
		
		output[[distance_decay_plot_ID]] <- renderPlot({
			isolate({
				plot(dist_decay_quadrats(sbssampling_quadrats(), method = "bray", binary = F))
			})
		})

      # make a note of the ID of this section, so that it is not repeated accidentally
      sbsvalues[[divID]] <- TRUE
      
      # create a listener on the newly-created button that will remove it from the app when clicked
      observeEvent(input[[btnID]], {
        removeUI(selector = paste0("#", divID))
        
        sbsvalues[[divID]] <- NULL
        
      }, ignoreInit = TRUE, once = TRUE)
      
      # otherwise, print a message to the console
    } else {
      message("The button has already been created!")
    }
  })

}) # end of server()





dist_decay_quadrats <- function(samples1, method = "bray", binary = F)
{
   com_mat <- samples1$spec_dat[rowSums(samples1$spec_dat) > 0,]
   d <- stats::dist(samples1$xy_dat[rowSums(samples1$spec_dat) > 0,])

   similarity <- 1 - vegan::vegdist(com_mat, method = method,
                                    binary = binary)
   similarity[!is.finite(similarity)] <- NA

   dat_out <- data.frame(distance = as.numeric(d),
                         similarity = as.numeric(similarity))

   # order by increasing distance
   dat_out <- dat_out[order(dat_out$distance), ]

   class(dat_out) <- c("dist_decay", "data.frame")

   return(dat_out)
}