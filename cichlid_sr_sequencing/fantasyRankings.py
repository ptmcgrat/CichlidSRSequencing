import pandas as pd

data = pd.read_html('https://www.nbcsports.com/nfl/matthew-berry/news/matthew-berrys-latest-positional-rankings-for-2024-fantasy-season', header = 0)
data2 = pd.read_html('https://www.nbcsports.com/nfl/matthew-berry/news/matthew-berrys-overall-top-200-for-2024-fantasy-football-season', header = 0)

output_file = 'FantasyRankings.xlsx'

with pd.ExcelWriter(output_file) as writer:  
	data[0].to_excel(writer, sheet_name = 'QB')
	data[1].to_excel(writer, sheet_name = 'RB')
	data[2].to_excel(writer, sheet_name = 'WR')
	data[3].to_excel(writer, sheet_name = 'TE')

	data2[0].to_excel(writer, sheet_name = 'Overall')
