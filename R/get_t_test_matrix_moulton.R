get_t_test_matrix_moulton <-
function(treatment_factor,data,cluster_factor,...)
{
	treatment_factor = as.factor(treatment_factor)
	l=levels(treatment_factor)
	t_test_matrix = matrix(nrow=length(l), ncol=length(l))
	dimnames(t_test_matrix)=list(as.character(l),as.character(l))
    moulton_factor_matrix = t_test_matrix
	for(a1 in l)
	{
		for(a2 in l)
		{
			if((length(data[treatment_factor==a1])==1) | (length(data[treatment_factor==a2])==1))
			{
				pval=1
                mf=NA
			}
			else
			{
				# t.test gives an error if all values are constant and equal 
				if(all(data[treatment_factor==a1] == (data[treatment_factor==a1][1])) & all(all(data[treatment_factor==a2] == (data[treatment_factor==a1][1]))))
				{
					pval=1
                    mf=NA
				}
				else
				{
					t_test=t.test.moulton(x=data[treatment_factor==a1],y=data[treatment_factor==a2],
                        cluster_x=cluster_factor[treatment_factor==a1],cluster_y=cluster_factor[treatment_factor==a2],...)
					pval=t_test$p.value
                    mf=t_test$moulton_factor
                    
				}
			}
			t_test_matrix[as.character(a1),as.character(a2)]=pval
            moulton_factor_matrix[as.character(a1),as.character(a2)]=mf
		}
	}
    
    attr(t_test_matrix,"moulton_factor")<-moulton_factor_matrix
	
	return(t_test_matrix)
	
			
	
}

