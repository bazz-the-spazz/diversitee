## Function to calculate diversities
		diversitee <- function(x, q, margin=1, effective=TRUE, na.rm=TRUE, se=F, allow.incorrect.beta=FALSE){
			# based on Jost 2006 & 2007 calculate the (effective) diversity of a community. Loosely adapted from the diversity function in vegan.

			# Organize data
			x <- drop(as.matrix(x))
			if (!is.numeric(x)) stop("input data must be numeric")
			if (any(x < 0, na.rm = TRUE)) stop("input data must be non-negative")

			# Function to Calculate Effective Diversity
			f <- function(x, Q, na=na.rm, raw=FALSE){

				# Get Proportional Data
				x <- x[x>0] # Zeros are ignored
				if(na) x <- x[!is.na(x)]
				p <- x/sum(x)

				if(NA %in% p){
					data <- NA
				} else {
					if(Q==1){
						data <- exp(-sum(p*log(p))) #q=1 is not definded, calculate Shannon
						if(raw) data <- -sum(p*log(p)) # needed for Alpha
					} else {
						data <- sum(p^Q)^(1/(1-Q)) # Else use Jost's universal formula
					}
				}
				return(data)
			}


			# Apply the before defined function f to the data(frame)
			if(length(dim(x))>1){
				D <- apply(x, MARGIN = margin, FUN = function(x)f(x, Q=q) )
			} else {
				D <- f(x,Q=q)
			}


			# From Jost 2007 :PARTITIONING DIVERSITY INTO INDEPENDENT ALPHA AND BETA COMPONENTS
			# True Alpha diversity
			if(length(dim(x))>1){
				## First Calculate weights
				w <- apply(x, MARGIN = margin, FUN = function(x)f(x, Q = 0))
				w <- w/sum(w) # w is the statistical weight of Community j (usually the number of individuals in Community j divided by the total number of individuals in the region)

				if(q!=1){
					if(length(unique(w))==1 | allow.incorrect.beta){ # if weights w are not equal only Shannon might be used!
						sums <- D^(1/(1/(1-q)))
						A <- ( sum((w^q * sums))/ sum(w^q) )^(1/(1-q)) # eq 11a
					} else {
						A <- NA
					}
				}
				if(q==1){ # use formula 11b for q=1 (Shannon)
					lamdas <- apply(x, MARGIN = margin, FUN = function(x)f(x, Q=q, raw=TRUE) )
					A <- exp(sum(w*lamdas))
				}

				# Gamma? the numbers equivalent of the diversity index of the pooled samples
				G <- f( t(colSums(x/rowSums(d, na.rm = T), na.rm = T))/nrow(x), Q=q)



				# Beta (Following from eq9)
				B <- G/A



				# Standard Error

				if(se & (q==1 | allow.incorrect.beta)){
					se.a <- se.b <- se.g <- numeric()
					for(i in 1:(nrow(x)-1)){
						for(j in (i+1):nrow(x)){
							xx <- x[c(i,j),]
							w <- apply(xx, MARGIN = margin, FUN = function(x)f(x, Q = 0))
							w <- w/sum(w)
							se.a <-  c(se.a, exp(sum(w*apply(xx, MARGIN = margin, FUN = function(x)f(x, Q=q, raw=TRUE) ))))
							se.g <- c(se.g, f( t(apply(xx, MARGIN = ifelse(margin==1,2,1) , sum, na.rm=T)) ,Q=q))
							se.b <- c(se.b, se.g[length(se.g)]/se.a[length(se.a)])
						}
					}
					se.a <- sd(se.a)/sqrt(ifelse(margin==1, nrow(x), ncol(x)))
					se.b <- sd(se.b)/sqrt(ifelse(margin==1, nrow(x), ncol(x)))
					se.g <- sd(se.g)/sqrt(ifelse(margin==1, nrow(x), ncol(x)))
				}



			}


			# Transform to non-effective Indices (if needed)
			if(effective==FALSE){
				if(q==0){
					D <- D
					I <- "Species Richness"
				}
				if(q==1){
					D <- log(D)
					I <- "Shannon entropy"
				}
				if(q==2){
					D <- 1-(1/D)
					I <- "Gini-Simpson index"
				}
				if(!(q %in% 0:2)) {
					Value <- NA
					I <- NA
				}
				return(list(Value=D, index=I, gamma=G, alpha=A, beta=B))
			} else {
				if(se){
					return(list(D=D, gamma=G, alpha=A, beta=B, se.g=se.g, se.a=se.a, se.b=se.b))
				} else {
					return(list(D=D, gamma=G, alpha=A, beta=B))
				}
			}
		}
