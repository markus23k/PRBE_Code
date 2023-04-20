wl2.df %>% filter(day == 60) %>% mutate(sum_iotas = iota * Window.length) %>% 
  select(4) %>% unlist() %>% c() %>% diff() %>% unname()


great_idea <- data.frame()

for(i in 1:60){
  great_idea <- rbind(great_idea, 
                      wl2.df %>% filter(day == i) %>% 
                        mutate(sum_iotas = iota * Window.length, day_est = day + Window.length,
                               iota_est = lead(sum_iotas) - sum_iotas))
}

plot(great_idea %>% group_by(day_est) %>% summarise(Mean = mean(iota_est, na.rm = T)))
 
