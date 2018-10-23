# running models to extract vital rate parameters for Early and Late phase
#### Survival ####+++++++++++++++++++++++++++++++++++++++++++++++++++++++
e.survival <- lme4::glmer(survival ~ 1 +
                            + (1|FishID) + (1|reach) + (1|sampling),
                          data=data.early,
                          family=binomial,
                          lme4::glmerControl(optimizer=c('bobyqa')))
# in early phase survival is not affectd by density
m.survival <- lme4::glmer(survival ~ s_loc_density +
                            + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                          data=data.mid,
                          family=binomial,
                          lme4::glmerControl(optimizer=c('bobyqa')))
# in mid survival is positively affected by density
l.survival <- lme4::glmer(survival ~ s_loc_density +
                            + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                          data=data.late,
                          family=binomial,
                          lme4::glmerControl(optimizer=c('bobyqa')))
# in late phase survival is positively affected by density

#### Recruitment ####+++++++++++++++++++++++++++++++++++++++++++++++++++++++
e.recruit <- glmmTMB::glmmTMB(recruited ~ s_loc_density +
                                + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                              ziformula = ~ 1 +
                                + (1|FishID) + (1|reach)  + (1|sampling),
                              family="truncated_poisson"(link = "log"),
                              data=data.early,
                              na.action = na.omit)
# in the early phase, density affects the number of offspring but not the probability to produce any.
m.recruit <- glmmTMB::glmmTMB(recruited ~ s_loc_density +
                                + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                              ziformula = ~ s_loc_density +
                                + (1|FishID) + (1|reach)  + (s_loc_density|sampling),
                              family="truncated_poisson"(link = "log"),
                              data=data.mid,
                              na.action = na.omit)
# in the mid phase both reproduction probability and number of offspring are negatively affected by density
l.recruit <- glmmTMB::glmmTMB(recruited ~ 1 +
                                + (1|FishID) + (1|reach) + (1|sampling),
                              ziformula = ~ 1 +
                                + (1|FishID) + (1|reach)  + (1|sampling),
                              family="truncated_poisson"(link = "log"),
                              data=data.late,
                              na.action = na.omit)
# the late phase is not affected by density

#### Adult dispersal probability ####++++++++++++++++++++++++++++++++++++++++++++++++
e.dispersal.p <- lme4::glmer(p_moving ~ 1 +
                               + (1|FishID) + (1|reach) + (1|sampling),
                             family = 'binomial',
                             data = data.early,
                             lme4::glmerControl(optimizer=c('bobyqa')))
# no effect of density on dispersal in early phase
m.dispersal.p <- lme4::glmer(p_moving ~ s_loc_density +
                               + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                             family = 'binomial',
                             data = data.mid,
                             lme4::glmerControl(optimizer=c('bobyqa')))
# positive density dependent dispersal probability for the mid part (higher dispersal at high density)
l.dispersal.p <- lme4::glmer(p_moving ~ s_loc_density +
                               + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                             family = 'binomial',
                             data = data.late,
                             lme4::glmerControl(optimizer=c('bobyqa')))
# strong positive effect of density on dispersal in the late phase

#### Adult dispersal distance ####+++++++++++++++++++++++++++++++++++++++++++++++++++++++
e.dispersal.dist <- lme4::lmer(log(dist_moved) ~ 1 +
                                 + (1|FishID) + (1|reach) + (1|sampling),
                               data = data.early)
# early phase distance does not depend on density
m.dispersal.dist <- lme4::lmer(log(dist_moved) ~ s_loc_density +
                                 + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                               data = data.mid)
# neither does the mid phase
l.dispersal.dist <- lme4::lmer(log(dist_moved) ~ 1 +
                                 + (1|FishID) + (1|reach) + (1|sampling),
                               data = data.late)
# and neither the late phase


#### Natal dispersal probability ####++++++++++++++++++++++++++++++++++++++++++++++++
e.natal.p <- lme4::glmer(p_natal ~ 1 +
                           + (1|Maternity) + (1|reach.birth) + (1|sampling.birth),
                         family = 'binomial',
                         data = data.early,
                         lme4::glmerControl(optimizer=c('bobyqa')))
# no density dependence early on
m.natal.p <- lme4::glmer(p_natal ~ s_loc_density.mid +
                           + (1|Maternity) + (1|reach.birth) + (s_loc_density.mid|sampling.birth),
                         family = 'binomial',
                         data = data.mid,
                         lme4::glmerControl(optimizer=c('bobyqa')))
# strong negative density dependence in the middle
l.natal.p <- lme4::glmer(p_natal ~ s_loc_density.mid +
                           + (1|Maternity) + (1|reach.birth) + (s_loc_density.mid|sampling.birth),
                         family = 'binomial',
                         data = data.late,
                         lme4::glmerControl(optimizer=c('bobyqa')))
# negative effect late 


#### Natal dispersal distance ####++++++++++++++++++++++++++++++++++++++++++++++++
e.natal.dist <- lme4::lmer(log(dist_natal) ~ 1 + 
                             + (1|Maternity) + (1|reach.birth) + (1|sampling.birth),
                           data = data.early)
# no effect of density on distance in early phase
m.natal.dist <- lme4::lmer(log(dist_natal) ~ s_loc_density.mid + 
                             + (1|Maternity) + (1|reach.birth) + (s_loc_density.mid|sampling.birth),
                           data = data.mid)
# higher distance at higher density in the mid phase
l.natal.dist <- lme4::lmer(log(dist_natal) ~ 1  + 
                             + (1|Maternity) + (1|reach.birth) + (1|sampling.birth),
                           data = data.late)
# no effect in late