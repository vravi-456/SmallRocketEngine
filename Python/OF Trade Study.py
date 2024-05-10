# OF = 2 is fuel rich compared to OF that maximizes Isp, more fuel available near walls to create cooler combustion zones
OF_list = np.linspace(0.1, 5, 100, True)
cstar_list = []
for OF in OF_list:
  c_star = obj.get_Cstar(P_c, OF)
  cstar_list.append(c_star)

plt.plot(OF_list, cstar_list)
i = np.argmax(cstar_list)
OF_list[i]