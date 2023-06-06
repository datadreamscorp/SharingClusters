### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 5cdc229c-48b8-11ec-04cf-8dec64279c6f
begin
	using PlutoUI
	
	md"""
	
	#### The evolution of sharing networks under conditions of risk and uncertainty
	
	##### A project-in-progress by Alejandro Pérez Velilla
	
	###### Motivation
	
	Welcome. This is a work-in-progress on the relationship between environmental risk, uncertainty, reciprocal sharing norms and social structure. Most accounts of human behavior, culture and group dynamics emphasize the ways in which groups working in cooperative and coordinated ways can achieve high payoffs to behavior and promote beneficial results like innovation. Another perspective that might be worth considering more closely is the way in which organizing in groups may aid risk mitigation when in changing environments and uncertain conditions (which is something all of us experience in life). Reciprocal sharing practices (known widely in the anthropological literature as "reciprocity"), in which individuals enter lasting relationships of reciprocal exchange of beneficial goods have been observed in societies around the world, and they are commonly accompanied by systems of norms that implicitly or explicitly mediate them. The following pamphlet shows the results of using a game-theoretic framework to explore the conditions in which reciprocal sharing practices can invade populations of non-sharers, and the effect that they might have on social structure once they have taken hold.
	
	There are two parts to this work. On one side there is the evolution of sharing practices, which will require careful examination of the ways in which we evaluate the cultural fitness of social strategies. On the other side we have the social consequences of sharing and the group dynamics it leads to once it has successfully invaded.
	
	I have included a reference guide for those who are curious about the mathematical formulation of the model that I used to derive the main results. I have kept it sparse because in this pamphlet I want to concentrate on a qualitative understanding of the dynamics rather than their mathematical derivation. For that reason, if you are not interested in the math, you may simply skip the next section. If you are, I apologize for the sparseness, if you want a more complete run-down please do contact me. Overall, there will be minimal math throughout this presentation of the work.
	
	###### Model reference guide (for the math-curious)
	
	 $B$: Benefits from successful work at niche or finding a resource patch. 
	
	 $u \in (0,1)$: Individual security. Individuals' rate of success of niche work / patch search.
	
	 $\mu \in (0,1)$: Group security. Additional security obtained for being part of the group. It may come from group-specific dynamics such as information-sharing.
	
	 $\rho \in (0,1)$: Total security. Individual security plus group security.
	
	 $s \in [0,1]$: Sharing coefficient. Proportion of earned benefits destined for sharing.
	
	 $k_i$: Degree of individual $i$.
	
	 $\bar{k}$: Expected degree of the group.
	
	 $C$: Per-connection cost of maintaining social interaction.
	
	 $z$: Minimum proportion of benefit for individual survival.
	
	 $N$: Number of individuals in the sharing cluster.
	
	 $\delta \in [0,\infty)$: risk preference. 
	
	**Strategy expected payoffs ("objective" framework)**:
	
	$V_i(Loner) = uB$
	
	$V_i(Informer) = \rho B - k_i C$
	
	$V_i(Sharer) = (1 - s) uB + s u k_i \frac{B}{\bar{k}} - k_i C$
	
	$V_i(Sharer-Informer) = (1 - s) \rho B + s \rho k_i \frac{B}{\bar{k}} - k_i C$
	
	**Strategy payoffs ("subjective" framework, perfect information)**:
	
	$V_i(Loner) = u^{\delta} B$
	
	$V_i(Informer) = \rho^{\delta} B - k_i C$
	
	$V_i(Sharer) = (1 - s) u^{\delta} B + s u k_i \frac{B}{\bar{k}} - k_i C$
	
	$V_i(Sharer-Informer) = (1 - s) \rho^{\delta} B + s \rho k_i \frac{B}{\bar{k}} - k_i C$
	
	###### The model
	
	Imagine that you are going through life looking for opportunities to obtain things that are beneficial to you, while in an environment that is constantly changing in ways that are not completely predictable to you (I know, crazy, right?). Every time period, you have to go out there and find yourself some niche (or patch, as they are commonly known in the foraging literature) and try to extract some benefit from it. Sometimes you find a good niche, and you are able to extract some benefit from it (for example, you apply for a job and get accepted). Some other times, you are not as lucky, and the search yields nothing (the CVs you sent out weren't even looked at). Jobs don't last forever, so every time period you have to get in the grindset and look for something new (or work at something you know, with uncertain returns).
	
	Your uncertain search for opportunities is characterized by an important aspect: risk. You don't know whether you'll get the job or not, although you may have some idea of what your chances are. If you know exactly what your chances are, then you know what your rate of success (which we will also call "security") is. Regardless, it is trivially true to say that whatever your chances are, you will always want them to be higher. We all want complete certainty, or whatever is closest to it, when life puts us through the daily hustle. And for that, maybe going at it alone is not the best strategy. Maybe grouping up with other people to inform each other of possible opportunities can better our chances, even if connecting with others entails some cost (because we have to spend time and resources to maintain our communications). I am better insured when I can count not only on my chances of finding a job opportunity on my own, but also on my cousin Roman calling me when there are extra job openings at the taxi company he just got hired at. And Roman enjoys the benefit of counting on me when I find a viable opportunity with enough space for us both.
	
	Let us imagine now that we are part of a population where everyone is going at it alone. Every time period people go around, look for opportunities, exploit them when they find them, and get nothing if they don't. We call this the loner strategy. But at some point, some of us decide to get together and start sharing information with one another, forming an information-sharing network. This is the informer strategy. In game-theoretic terms, assuming (for simplicity) that I'll always get a job that I'm informed about, getting together in networks of interaction and informing one another of opportunities will be more favorable than going at it alone if the following condition is true:
	
	$(1 - u) B > \bar{k} C$
	
	This expression says that as long as I am not making too many connections, or these are not too costly (with respect to benefits) to make in the first place, then the higher the risk, the better it will be to get together and share information. Thus, informers will tend to get enough connections to balance out the increased security they get from being part of an information sharing network and the costs of connecting to other individuals.
	
	"""
end

# ╔═╡ 8514d31b-b9e9-4e45-bb65-28b19d83d637
LocalResource("loners_informers.png")

# ╔═╡ 45e02182-d3ff-4ad8-8444-476fc7c4f14e
md"""

So if a group of informers gets together and starts doing their thing, they will be enhancing their chances of getting opportunities to exploit in this uncertain environment, and the higher the risk of not getting anything out of our search, the more advantageous this strategy will be. Surely others will quickly notice that some people are doing better by adopting this strategy, and will choose to adopt them themselves. In cultural evolutionary terms, the informer strategy has a higher cultural fitness than the loner strategy; social learners in the population will find it more attractive than the loner strategy when making the decision of what strategy to choose or change to.
	
What about reciprocal sharing practices? If, when I get into a social relationship with another person, we agree to not only share information but also share a part of our earned benefits with one another, will that mean I will be better insured? To find out, we do the math. We call the strategy that gets together in networks only to share benefits (but not information) the sharer strategy. And the sharers that do share information as well are using the sharer-informer strategy.

"""

# ╔═╡ 493f7e64-a6ec-4b13-b990-f7cc7a00550b
LocalResource("sharers_lose.png")

# ╔═╡ f381c649-692c-47df-9bc5-97960e3d3526
md"""

Doing the calculations, it turns out that sharers are at a disadvantage even against loners (because they pay a cost to connect with one another, but their chances of getting benefits from niche work are the same as loners') and sharer-informers have the exact same payoff as plain informers. Why would we even bother sharing with each other if it is not going to do anything for us?
	
Why would reciprocal sharing ever take hold in conditions like these?
	
###### Changing perspectives: from all-knowing gods to boundedly-rational mortals
	
The normal way in which we have been calculating the cultural fitness of strategies is analogous to what in Economics is known as Expected Utility Theory (EUT). This is also popular practice in decision theory and risk management, but it is not without critics. The problem with this formulation of decision problems is that it assumes that individuals who make choices by evaluating their possible outcomes are not sensitive to risk factors like payoff volatility in time (also known as a strategy's variance). Instead, most of the decision-making weight in this framework assumes we calculate how well we'll be doing *in the long run*, and we don't really care for what happens in any particular time period. And as we saw above, in the long run, sharers get even less benefits that loners (because they're paying for connections when loners aren't) while sharer-informers have the exact same expected payoff of an equivalent cluster of informers. But this turns out to not be true if we do care about volatility.

"""

# ╔═╡ 9efecfe8-2890-4711-989b-7ed1684a7d2e
LocalResource("risk_aversion.png")

# ╔═╡ 27e2f4ba-1b4a-412e-87ae-e55408b70a66
md"""

While a full treatment of why we might want to avoid payoff volatility is beyond this pamphlet, an intuitive way that we might think about this subject is by asking ourselves: when will I prefer strategies that will provide me with moderate but secure payoffs in time, versus strategies that might provide very high payoffs sometimes but very low payoffs other times? There are many situations in which we want the peace of mind that comes with low volatility. For example, paying for life investments like education might require us to have a steady source of income. An even higher-stakes scenario is childcare: a steady source of income is necessary to provide for a child, especially during critical stages of development, which may have a strong impact on the child's future. The point is that for many human affairs, a steady moderate payoff will not be the same, from a decision-making perspective, as the possibility of getting a high payoff accompanied by the possibility of getting a low payoff (or, in this simplified model, nothing). By being blind to this distinction, the standard mathematical framework does not allow for the evolution of sharing. We require a change in perspective. For this, we have to consider the state of knowledge a social learner is in at the moment of deciding which target individual to learn from, as well as the learner's preferences regarding payoff volatility in time. In modeling terms, we have to think about how the observed payoffs and the information we have about them translate into cultural fitness. We need to move from an objetive point-of-view to a subjective one.

"""

# ╔═╡ 55a8e78e-fb80-4758-be25-4c613b2a96d7
LocalResource("change_of_perspective.png")

# ╔═╡ db56480f-14fa-4403-beb9-2325160a5f3f
md"""

We will focus on the risk preference for what remains of this work, so we will assume people's estimates about the rate of success of strategies is accurate (e.g. they look like the actual, objective rates of success). Nevertheless, I believe rate of success estimation and its consequences are very rich subjects that could have important consequences for general social learning dynamics.

###### The evolution of sharing and the group-level consequences of risk aversion

Under this new framework, in which the social learner's risk preference plays a part, sharing can now invade. The condition for the invasion of the sharing-informing strategy in a group of informers is

$s \rho (1 - \rho^{\delta - 1}) B > 0$

When this condition is satisfied, a networked cluster of sharer-informers will do better, in the eyes of social learners, than plain informers. As a consequence, learners will choose to join the network where sharing is also happening instead of the one where individuals are only sharing information. Similar conditions hold for the invasion of simple sharing and sharing-informing over the loner strategy. Interestingly, this condition requires that individuals are averse to payoff volatility ($\delta > 1$, also known as risk-averse or "pessimistic"). When risk-averse, the possibility of low payoffs negatively affects our perception of higher-risk strategies, even when there's a chance we'll get a high payoff in a single time period.

"""

# ╔═╡ 06891540-3486-4ccc-850d-d8d7b09d61b5
LocalResource("sharing_variance.png")

# ╔═╡ dbb03577-87ee-41a1-b5ec-3e6321b94e40
LocalResource("trajectories.png")

# ╔═╡ 3025f2a1-2b13-42f4-a494-043df57384cf
md"""

Another way of looking at the problem is in the case where stakes are particularly high, and getting a non-zero minimum amount of benefit is a true necessity, then individuals will be maximally risk-averse and the loner and informer strategies become totally undesirable. Meanwhile, the sharing strategies can not only face these conditions, but it can be shown that they can lead to group survival.

"""

# ╔═╡ fd1be83c-9102-48b9-aae0-23e65652f82c
LocalResource("group_survival.png")

# ╔═╡ a8081ba4-360d-4a0e-9867-392baa8e5d19
md"""

The nicest part of this result is that risk aversion is a common subject in the economic and decision-theoretic literature, due to it being repeatedly observed in empirical settings. Peoples' relative aversion to risk has been the subject of much work and discussion, with arguments attributing its common presence in our psychologies to proposed genetic factors arising out of adaptationist hypotheses, as well as to developmental, demographic and cultural factors related to social learning and environmental differences (Malek & Einsenhauer, 2001). In the evolutionary literature, risk aversion/pessimism has recently been the subject of theoretical inquiry, showing that in the case of biological fitness in age-structured populations, employing pessimistic decision-making can lead to fitness maximization (Price & Jones, 2020). Moreover, more complex modeling work in the anthropological literature on human foraging provide further support for the link between risk minimization and sharing practices in stochastic environments (Winterhalder, 1986; Winterhalder, Lu & Tucker, 1999).

The model we're examining here has a clear conceptual link with the above-mentioned literature, while examining the dynamics of sharing in a cultural evolutionary scenario where the decision-makers are the social learners. It also has the benefits of parsimony, mathematical tractability and easy interpretability.

So reciprocal sharing can evolve in a cultural population. That is nice, but is it the full story? The answer is **not at all**. Once sharing has invaded, sharing networks will inevitably grow with respect to loners and informer networks. Where does this growth take us? Turns out the answer is far from straightforward, and it reveals a side of sharing that is not cooperative, but competitive. This will push us into the territory of cultural anthropology, where many reciprocity-fueled social systems around the world have been given careful observation and analysis.

"""

# ╔═╡ 0b967efc-7bdc-4de6-98c6-14aafa23e1e9
md"""

###### Sharing networks' internal competition, group destabilization and fission

If we are a social learner that has decided to imitate a sharer's strategy, we are adquiring a bundle of cultural attributes. We are accepting to conform to a sharing norm with our network peers, and to provide them with information on available exploitable niches. We are quite possibly adopting social markers that identify us as part of the sharing group and that let us coordinate ourselves with other sharers we might meet. But, crucially, we will likely also identify the target individual's sociability as a non-trivial part of their strategy. Why's that? Because it turns out that the number of network peers a sharer accumulates has an impact on the shared benefits they receive. The more peers one has, the more shared benefits one can accumulate. This means that, within the sharing cluster, we are more likely to imitate high-degree individuals, and so we will go ahead and make around the same number of connections that they exhibit. In other words, we will imitate their sociality patterns, beyond just the norms and markers necessary for coordination.

This sets off a counter-intuitive feedback loop. As new sharers enter the network, they will seek to imitate the sociality patterns of high-degree agents, and so they will make connections that mirror the ones of agents in the tail-end of the network's degree distribution. This inevitably leads to a more connected network. More connections, more sharing, right? What's not to love? The problem with this is clear when we use a toy example.

"""

# ╔═╡ 054c039e-3e16-472f-a80d-9c4b460660a8
LocalResource("internal_competition.png")

# ╔═╡ 6107aaaf-ee76-4a01-9f1c-ebc1220d96b4
md"""

Let's say I belong to a big club with many members who don't necessarily all know each other, and in which I have made friends with whom I periodically meet one-on-one and share quality time with. We can see these relationships as social pacts in which I reciprocally share free time with friends. But if my friends suddenly start getting new friends (perhaps because the club is popular and it attracts new members who start befriending previous club members), the share of time they can dedicate to me one-on-one becomes less and less, because free time is limited. I am putting time aside for friends that now don't have much time for me, which is kind of depressing, so I choose to imitate the club members who seem to be doing best and look for more friends in the club to make up for the loss. Let's say I do so. Since everyone is trying to split their one-on-one time with their friends equally, my new friends' previous friends will inevitably perceive a loss of shared time with them, which will push them to get more friendships in the club to make up for it. In other words, low network degree individuals can be quite inconvenienced by high degree individuals accumulating ties. This drives them to imitate high degree individuals and become high degree themselves. If left unchecked, this situation leads the club's expected degree to increase, which means that a person entering the club seeking to form connections (or a current club member choosing to obtain more connections, which is formally equivalent in this framework) will increasingly perceive that the other club members that they can potentially befriend already have many friendships themselves, and thus will have very little time to dedicate to new friends. Because connections have a non-trivial cost of maintenance (perhaps because I live in a low-density suburban sprawl with little options of public transportation and high gas prices), there will be a point in this dynamic in which the mean degree of the network grows so large that the costs (both the current costs being payed and the new costs that getting more connections would entail) will overtake any benefits that I could be obtaining from this whole ordeal. At this point, I (and probably many others, with the possible exception of the highest-degree individuals) will seriously consider just leaving the club behind and concentrating on myself. Eventually, if enough people choose to abandon the club, the club's mean degree will decrease to a point in which new potential members might actually start seeing benefits to joining once more. The dynamic then repeats itself, leading the club's mean degree to cycle chaotically around a critical point: the point at which joining the club goes from beneficial to detrimental. In ecological terms, this critical point is something akin to a "population carrying capacity" but it applies to mean degree instead of population size.

"""

# ╔═╡ 765d8439-a3d7-4a8d-8939-b1850920dcde
LocalResource("crowding.png")

# ╔═╡ 98d1dda4-01ff-428b-b56f-3efb0d7b113c
md"""

The key is to realize that the network's degree is not directly controllable by individuals. In other words, there is a very real tension between an individual's want of reciprocity ties and the social consequences of forming them. Within the group, the individuals that acquire the most ties are at an advantage at the expense of individuals with comparatively less ties. At the group level, the competitive dynamics this asymmetry incentivizes lead to a drop in the average payoff of the group. In other words, getting greedy with the accumulation of reciprocity ties can seem like a good strategy, but it can eventually lead to an unstable group environment in which most people are not getting that much out of being part of the group, and thus are always on the verge of abandoning it.

What is even sadder is that this does not necessarily depend on greedy individuals accumulating all the reciprocity ties they can: if the club is constantly growing and people are making friendships even without trying to accumulate as much of them as they can, new group members have to make less than half the mean degree of connections on average for the mean degree of the club  to not increase in time. So even in situations where individuals are moderate with their tie formation our sharing dynamics are still likely to lead to group destabilization, just from the fact that new people are entering the group (this further shows the link between these dynamics and general population dynamics constrained by a carrying capacity).

But let's say I don't wait until things are at their worst to make a decision about leaving. After all, I am not the only one going through this: I talk to some of my closer friends, and we agree that things are getting a bit out of hand in this club. Maybe, just maybe, we should just get out now and start our own club, which will be less crowded and more chill.

As a way to reflect on how this simple abstract dynamic can be representative, at least qualitatively, of some of the group strategies and dynamics we may observe in real life, I will leave here two passages from independent ethnographic accounts of cooperation and village fissioning in Yanomamö society, which is, by all accounts, strongly organized around explicit reciprocity relations:

"Conflicts are reduced to a minimum within a lineage or co-resident group of agnates who stand, after all, as a virtual lineage (...) The phrase iba dibi sai (my real people) is an expression of this solidarity, as is the sharing of food and services, and the proximity of their dwellings and garden plots. The closeness and cooperation between co-resident agnates, made evident with a simple glance at a village's 'economic map,' are even more apparent when compared with the situation of people who have no agnates nearby. Agnatically isolated men may not be economically deprived, but they lack the comfort and security of counting on their relatives for whatever needs arise. Men who have no close agnates or live apart from them are more prone to accusations of theft, adultery, and other abusive acts. The phrase iba ai dibi (my others) denotes the ambivalence that exists in the relationship between affines and distant kin. They are both 'my people' (iba dibi) and 'other people'. The discomfort that accompanies much of their interactions can be observed in certain situations." (Ramos, 1995, p. 90)

"As villages grow larger, internal order and cooperation become more difficult, and eventually factions develop: Certain kin take sides with each other, and social life becomes strained. There appears to be an upper limit to the size of a group that can be cooperatively organized by the principles of kinship, descent, and marriage (...) kinship-organized groups can only get so large before they begin falling apart—fissioning into smaller groups." (Chagnon, 2012, p. 78)

In these passages, we can see that Yanomamö village fission can be brought on by demographic growth: as villages grow larger, kinship ties (the strongest basis of reciprocity) grow more diffuse, and individuals choose to prioritize their closest kinship ties in key interactions, which include reciprocal sharing. This leads to a degradation of village "solidarity" (to use Chagnon's term) that leads to within-village conflict and, eventually, village fission along close kin lines.

###### Hierarchy, stratification and the emergence of social classes

In 1963, the cultural anthropologist Marshall Sahlins authored a paper (today considered a classic of the anthropological literature) in which he examines precisely the ways in which reciprocity may lead to competitive dynamics such as the ones described above, by comparing the ways in which normative reciprocal sharing practices structure social and political life for societies in Melanesia and Polynesia (Sahlins, 1963). In his examination of the buildup of political power common in Melanesian societies, he stressed the role of reciprocity in the emergence of status hierarchies headed by Big-Men: leaders who have gained status through the accumulation of reciprocity ties inside their groups, and that employ them strategically in a form of competitive reciprocity between groups. You can imagine where this is going.

"One side of the Melanesian contradiction is the initial economic reciprocity
between  a  center-man  and  his  followers. For  his  help  they  give  their  help,
and  for  goods  going  out  through  his  hands  other  goods  (often  from  outside
factions)  flow  back  to his followers  by  the  same  path. The  other  side  is that a  cumulative  build-up  of  renown  forces  center-men  into  economic  extortion of  the  faction. Here  it  is  important  that  not  merely  his  own  status,  but  the standing  and  perhaps  the  military  security  of  his  people  depend  on  the  big-man's  achievements  in public  distribution. Established  at the head  of  a sizeable  faction,  a  center-man  comes  under  increasing  pressure  to  extract  goods from  his followers, to delay reciprocities owing them,  and to  deflect  incoming goods  back  into  external  circulation. Success  in  competition  with  other  big-men  particularly  undermines  internal-factional  reciprocities:  such  success  is precisely  measurable  by the  ability  to  give outsiders  more   than  they  can possibly  reciprocate. In  well  delineated  big-man  polities,  we  find  leaders  negating  the  reciprocal  obligations  upon  which  their  following  had  been  predicated.  Substituting extraction for  reciprocity, they must compel their people to  "eat  the  leader's  renown,"  as  one  Solomon  Island  group  puts  it,  in  return for  productive  efforts. Some  center-men  appear  more  able  than  others  to dam  the  inevitable  tide  of  discontent  that  mounts  within  their  factions,  perhaps  because  of  charismatic  personalities,  perhaps  because  of  the  particular social  organizations  in  which  they  operate.15 But  paradoxically  the  ultimate defense  of the center-man's  position  is some slackening  of his drive to enlarge the  funds  of  power. The  alternative  is  much  worse. In  the  anthropological record  there  are  not  merely  instances  of  big-man  chicanery  and  of  material deprivation  of  the  faction  in  the  interests  of  renown,  but  some  also  of  overloading  of social  relations  with  followers:  the  generation  of antagonisms, defections,  and  in  extreme  cases  the  violent  liquidation  of  the  center-man." (Sahlins, 1963)

While modeling the explicit dynamics of Big-Man reciprocity and status competition would probably require a slight alteration of the current model, the conclusions are the same, at least in a qualitative sense: the driving-up of competitive reciprocity leads to dynamics that are fundamentally unstable when not kept in check by other forces. The emergent hierarchy of internal competition in reciprocity networks will not lead to stable hierarchical patterns in time. Rather, it can promote, same as demographic growth in the above examples from Yanomamö society, an eventual destabilization of group solidarity.

So, going back to our model, let's say that people are a bit choosier with potential reciprocity peers. They will not want to connect with someone who will give them a shared benefit that they perceive as too low. Not only that, but they assume other people are choosy too. Why invest in someone who will not invest significantly in me? Better to concentrate on the ties I have that do share with me in a way I think is significant enough. Especially considering that with each additional tie I make, I thin out the shared benefits I confer to my other ties. In other words, I can only make a limited number of ties before people will start choosing to break up with me. So I better choose well.

In a scenario such as this one, we need only introduce a little heterogeneity in niches to obtain a different kind of hierarchy, a more stable one fueled by an uneven access to defensible resource niches with differing benefits and securities. Imagine now that the niches around the environment are not uniform: there are many "low-yield" niches that provide a decent amount of benefit with some risk, and a smaller amount of rare "high yield" niches that provide high amount of benefit with comparatively lower risk. Furthermore, assume that these niches are defensible, so that individuals that come to occupy them can essentially monopolize their use and regulate access to them.

"""

# ╔═╡ a42492e2-ff72-440f-9676-4091de573e76
LocalResource("niche_ownership.png")

# ╔═╡ 6916b493-9850-4c5c-8e28-8d5019b0bfe3
md"""

Since now individuals are mindful of the connections they are making, the model shows that it will always be better for an individual to prioritize connecting with the individuals in the high-yield niches first, as the shared benefits obtained by them are more secure (and maybe higher, but this need not be the case). Linking with these individuals, in sum, is a better personal insurance strategy than linking with individuals at low-yield niches, and so the higher the stakes are, the more individuals will prioritize connecting with high-yield individuals.

"""

# ╔═╡ 4d9a5210-eac4-45df-82c8-2f6ca078a7ee
LocalResource("tie_accumulation.png")

# ╔═╡ b526d505-4be4-440b-aeee-8ccb67939d52
md"""

High-yield individuals, on the other hand, have an additional upper hand in that they can invest in more connections than low-yield individuals, since the higher, more secure benefits that their niches provide them allows them to put aside more benefits for sharing. There is a limit, yes: this is only true as long as a high-yield individual does not accumulate too many connections that even their higher shared benefits become too diluted through their ego networks. However, this limit can be much higher than that of a low-yield individual. Mathematically, and only considering the case of high and low securities (which assumes equal potential benefits across all niches):

$k_j < \frac{u_H}{u_L} \bar{k}_L$

where $u_H$ is the security of high-yield individuals, $u_L$ is the security of low-yield individuals, $k_j$ is the degree of the focal high-yield individual $j$, and $\bar{k}_L$ is the mean degree of low-yield individuals. The higher the difference between the niche type risks, and the higher the mean degree of low-yield individuals is, the higher the number of connections a high-yield individual can comfortably make.

"""

# ╔═╡ 15c6e1e9-bff4-40dc-b067-c6837042d5e8
LocalResource("preferential_attachment01.png")

# ╔═╡ 0f434ab1-35a7-409b-be34-f8e38a4e7438
md"""

These hierarchies, backed up by defensible positions of environmental privilege, can be seen as a more stable alternative to the emergent competitive hierarchies that we saw above. This environmental privilege might come from actually occupying ecological niches, but also from enjoying a privileged position in societal niches. When surveying Polynesian chiefdoms, Sahlins writes:

"The  growth  of  a  political  system  such  as  the  Polynesian  constitutes
advance  over  Melanesian  orders  of  interpersonal  dominance  in  the  human
control  of  human  affairs. Especially significant for society  at  large  were
privileges  accorded  Polynesian  chiefs  which  made  them  greater  architects of funds  of  power  than  ever was  any  Melanesian  big-man. Masters  of their people and "owners" in a titular  sense of  group resources, Polynesian  chiefs  had  rights  of  call  upon  the  labor  and  agricultural  produce of  households  within  their  domains. Economic  mobilization  did  not  depend on,  as  it  necessarily  had  for  Melanesian  big-men,  the  de  novo creation  by the  leader  of  personal  loyalties  and  economic  obligations. A  chief  need  not stoop to  obligate  this  man  or  that  man,  need  not  by  a  series  of  individual  acts of generosity induce others  to support  him,  for economic  leverage  over a group  was  the  inherent  chiefly  due."

He continues:

"Redistribution  of  the fund  of power was the supreme art of Polynesian  politics. By  well-planned  *noblesse  oblige*  the  large  domain  of  a  paramount  chief  was held  together,  organized  at  times  for  massive  projects,  protected  against  other chiefdoms,  even  further enriched. Uses  of  the  chiefly  fund included  lavish hospitality and entertainments for outside chiefs and for the chief's own people,  **and succor  of  individuals  or  the  underlying  population at  large  in times  of  scarcities  —  bread  and  circuses**." (Sahlins, 1963. Emphasis is mine.)

"""

# ╔═╡ af1655b1-fab9-41fd-8a7f-8e7cd0bfa507
LocalResource("reciprocity_hierarchy.png")

# ╔═╡ 22cec7fe-0e33-4ee4-84c0-11cef42f6a02
md"""

Finally, it is worth giving a brief mention to what might happen in a situation in which there are several high-yield individuals in the population. Because of the obvious advantage in insurance, they will also prioritize connecting with other high-yield individuals before anyone else. If these elites keep growing in number and keep monopolizing high-yield niches, then at some point they will be making more connections with other high-yielders than with low-yielders. In this situation, only a few lucky low-yielders might be the ones that get to connect with high-yielders, who are mostly connecting with one another, and who only connect with low-yielders if they happen to have additional room for more connections once they have connected with each other. Facing this bleak panorama, most low-yielders will connect with other low-yielders, as that is still better than facing life alone. The result of this dynamic will be two well-differentiated clusters, only sparsely connected with each other, if at all. This is the emergence of social classes through reciprocity and risk-insurance. 

"""

# ╔═╡ 89977002-b834-48f7-8302-da6821ee11a8
LocalResource("preferential_attachment02.png")

# ╔═╡ 7429c089-8d88-4e89-8af7-1cac05ed335b
md"""

###### The possibility of egalitarianism

All this discussion of hierarchy out of reciprocity makes it seem like more egalitarian outcomes are not possible. But societies with strong egalitarian aspects do exist around the world, and they normally have an important element of reciprocity to them. Luckily, this modeling framework is flexible enough to represent egalitarianism in plausible terms: we only need to assume that individuals do not vary in their number of reciprocity ties. Instead, they all have the same amount of ties, because they are all connected to each other. In formal terms, the sharing network in this case becomes a fully-connected network, where everyone shares with everyone else.


"""

# ╔═╡ b514fa51-4b7f-480f-ac8f-71e42cf008c7
LocalResource("egalitarianism.png")

# ╔═╡ d8d92032-0f8b-45bf-9016-8f6149b9f84e
md"""

In this case, any possibility of hierarchy is completely flattened, because an individual's continued group membership is conditional on them accepting to share with all other group members. In exchange, the individual will also be shared with. No-one can gain more connections than anybody else, so no-one can get ahead of anybody else. The mathematical condition for this to succeed in a context of group survival is

$\bigg \{ (1 - z)\rho - z \bigg \} \frac{B}{N-1} > C$

Essentially, this says that an egalitarian group can achieve group survival if its productivity is high with respect to individuals' needs. Specifically, the benefits each individual potentially gets from niche work have to be at least more than double than what they need to survive. This is because only then there will be the chances of there being enough surplus so that individuals who did not get any benefit in the present time period can add enough shared contributions from successful group members in order to cover their needs. These members might end up earning benefits in the next time period, and they may end up ensuring the survival of individuals that aided their own survival in the previous time period. Also, the total security, that comes from individual and group securities, has to be high enough, which allows us to summarize the condition as: group productivity (magnitude of benefits and securities) must be high with respect to needs.

But the most crucial aspect is that the dynamics now truly depend on group size rather than mean degree: if the group grows too large, the benefits might once again thin out so much that people's obtained shared benefits might start getting dangerously low for individual (and, by extension, the group's) survival in time. Group size, on the other hand, is much easier to regulate than each person's individual number of network ties, a job which might be assumed by leaderships (like elders) and other political institutions. This dependence on group size in the model hints at two plausible insights.

First, that the higher productivity required by these egalitarian groups to sustain higher group size is likely to place them under strong constraints on scale. After all, keeping a group small in size not only ensures shared benefits do not get too thinned out, but also makes it easier for its members to gauge the growth of the group while taking into account what they are managing to get from their productive activities. If the group is growing, but productivity (again, in terms of security and niche benefits) isn't, at some point we are going to hit a limit of how many people we can have in our band without putting everyone else in danger. In this way, group decisions regarding group size and actions to be taken about it (like, for example, the decision to peacefully fission a large group into two smaller ones that occupy different environments as a remedy to too much growth) can be plausibly achieved with greater ease for smaller groups, as it is easier for group members (and leaderships) to keep tabs on everyone else. This might be a reason behind the observed common patterns of organization of egalitarian hunter-gatherer bands, in which individuals tend to be part of mobile bands of at most several dozen individuals (whose costs of connection are also negligible, because they occupy the same community space). However, there is nothing in the model denying the possibility of large-scale egalitarian societies.

Second, that egalitarianism at the level of benefits/utility might have a fitting companion to it in the form of political hierarchy and authority. This is because it is imaginable that the control over decisions regarding group size as well as an effective coordination of benefit redistribution throughout the group might require, at least in some cases, that an instituted leadership evaluates the group's conditions and takes quick decisions in favor of group welfare. The faster decisions have to be made (which might be a concern for groups in very unstable environments or subjected to very high stakes), the likelier it is that these groups will have to rely on some form of political authority to consolidate a fast, effective decision-making regarding potentially tough subjects. 

In any case, this modeling exercise has shown us that interdependence and group survival are achievable in egalitarian groups through risk-insurance, and that reciprocal sharing is a strategy through which this insurance might be carried out. I close this sections with another quote from Sahlins which I find relevant to the current discussion:

"On a very general view, the array of economic transactions in the ethnographic record may be resolved into two types. First, those "vice-versa" movements between two parties known familiarly as 'reciprocity.' The second, centralized movements: collection from members of a group, often under one hand, and redivision with this group. This is 'pooling' or 'redistribution'. On an even more general view, the two types merge. For pooling is an organization of reciprocities, a system of reciprocities - a fact of central bearing on the genesis of large scale redistribution under chiefly aegis." (Sahlins, 1972, p. 188)

I am glad to show that, even for a simple model like the one we have seen here, this conceptual merge does happen, and redistribution can indeed be seen as an organization of reciprocities, represented by a groups' network structure and the individual and societal constraints placed on it.

###### Discussion: ecological rationality and group dynamics

I want to close this pamphlet with a reflection on how considering how people learn from one another is, in its essence, a decision-theoretic problem. Incorporating this perspective into our models of cultural evolution will not only help us derive novel results at the group level, but it will help us tie these results to aspects of our individual psychologies and the factors that shape them.

As soon as we consider pessimism as an influential factor in the mix, reciprocal sharing is able to evolve. And as soon as it evolves, it exhibits qualitative group dynamics that match aspects of the diversity behind the grouping strategies (and the structures they lead to) that anthropologists and historians have extensively documented around the world, and some of the important theoretical considerations that they have sparked.

Human societies can be seen as massive risk-pooling systems. Under this perspective, strategies that seem "irrational" or unable to evolve under a standard point of view might instead prove to be effective and successful. Understanding how psychological factors like risk aversion (and the environmental factors that shape them throughout an individual's life) can serve as ways to frame proposed social learning strategies in terms of key concepts in decision theory. This can also help us move in a direction where we can unite these proposed social learning strategies (like homophily-biased or prestige biased learning) under a common framework of bounded rationality, instead of taking them for granted as evolved cognitive biases. By using pessimism as a filter for focusing on strategies that are not only observably well-paying, but also seem like they will be well-paying *for me*, then it might turn out that my choices are less because of some evolved bias and more because of my uncertainty about the environment and my use of social information (like social identity) in order to mitigate perceived risks. In this view, human social learning is seen as ecologically rational decision-making.

"""

# ╔═╡ a267b9c7-5108-4be6-816c-fd3bdd0d2417
LocalResource("conclusion.png")

# ╔═╡ 33f357b7-fdae-4ceb-be57-9217b8fdc497
md"""

###### References

Chagnon, N. A. (2012). The Yanomamo. Cengage Learning.

Halek, M., & Eisenhauer, J. G. (2001). Demography of risk aversion. Journal of Risk and Insurance, 1-24.

Price, M. H., & Jones, J. H. (2020). Fitness-maximizers employ pessimistic probability weighting for decisions under risk. Evolutionary Human Sciences, 2.

Ramos, A. R. (1995). Sanumá memories: Yanomami ethnography in times of crisis. Madison: The University of Wisconsin Press.

Winterhalder, B. (1986). Diet choice, risk, and food sharing in a stochastic environment. Journal of anthropological archaeology, 5(4), 369-392.

Winterhalder, B., Lu, F., & Tucker, B. (1999). Risk-senstive adaptive tactics: models and evidence from subsistence studies in biology and anthropology. Journal of Archaeological Research, 7(4), 301-348.

Sahlins, M. D. (1963). Poor man, rich man, big-man, chief: political types in Melanesia and Polynesia. Comparative studies in society and history, 5(3), 285-303.

Sahlins, M. D. (1972). Stone Age Economics (No. 306.3 S2).

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.19"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0bc60e3006ad95b4bb7497698dd7c6d649b9bc06"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "e071adf21e165ea0d904b595544a8e514c8bb42c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.19"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─5cdc229c-48b8-11ec-04cf-8dec64279c6f
# ╟─8514d31b-b9e9-4e45-bb65-28b19d83d637
# ╟─45e02182-d3ff-4ad8-8444-476fc7c4f14e
# ╟─493f7e64-a6ec-4b13-b990-f7cc7a00550b
# ╟─f381c649-692c-47df-9bc5-97960e3d3526
# ╟─9efecfe8-2890-4711-989b-7ed1684a7d2e
# ╟─27e2f4ba-1b4a-412e-87ae-e55408b70a66
# ╟─55a8e78e-fb80-4758-be25-4c613b2a96d7
# ╟─db56480f-14fa-4403-beb9-2325160a5f3f
# ╟─06891540-3486-4ccc-850d-d8d7b09d61b5
# ╟─dbb03577-87ee-41a1-b5ec-3e6321b94e40
# ╟─3025f2a1-2b13-42f4-a494-043df57384cf
# ╟─fd1be83c-9102-48b9-aae0-23e65652f82c
# ╟─a8081ba4-360d-4a0e-9867-392baa8e5d19
# ╟─0b967efc-7bdc-4de6-98c6-14aafa23e1e9
# ╟─054c039e-3e16-472f-a80d-9c4b460660a8
# ╟─6107aaaf-ee76-4a01-9f1c-ebc1220d96b4
# ╟─765d8439-a3d7-4a8d-8939-b1850920dcde
# ╟─98d1dda4-01ff-428b-b56f-3efb0d7b113c
# ╟─a42492e2-ff72-440f-9676-4091de573e76
# ╟─6916b493-9850-4c5c-8e28-8d5019b0bfe3
# ╟─4d9a5210-eac4-45df-82c8-2f6ca078a7ee
# ╟─b526d505-4be4-440b-aeee-8ccb67939d52
# ╟─15c6e1e9-bff4-40dc-b067-c6837042d5e8
# ╟─0f434ab1-35a7-409b-be34-f8e38a4e7438
# ╟─af1655b1-fab9-41fd-8a7f-8e7cd0bfa507
# ╟─22cec7fe-0e33-4ee4-84c0-11cef42f6a02
# ╟─89977002-b834-48f7-8302-da6821ee11a8
# ╟─7429c089-8d88-4e89-8af7-1cac05ed335b
# ╟─b514fa51-4b7f-480f-ac8f-71e42cf008c7
# ╟─d8d92032-0f8b-45bf-9016-8f6149b9f84e
# ╟─a267b9c7-5108-4be6-816c-fd3bdd0d2417
# ╟─33f357b7-fdae-4ceb-be57-9217b8fdc497
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
