# Multi-UAV Energy Consumption Simulator
The Multi-UAV Energy Consumption (MUAVE) Simulator is a Python-based simulator that allows to compute the energy consumption of multiple Unmanned Aerial Vehicles (UAVs), acting as Flying Access Points (FAPs) to provide wireless connectivity to Ground Users (GUs).

As input, it receives a .txt file named "GUs.txt" that identifies the number of groups of GUs, the total number of GUs, their spatial positions, and their traffic demand.

As output, it provides the energy consumption of the FAPs following 2 state-of-the-art algorithms and an energy consumption reduction percentual value when compared to a hovering state.
