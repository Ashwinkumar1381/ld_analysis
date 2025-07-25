#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float dt = 5e-4;
	int frameW = int(1e3);

	Trajectory *TRAJ = new Trajectory(dt ,frameW, "xyz");
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data29/traj2.xyz");
	
	TRAJ -> openTrajectory(true);

	TRAJ -> closeTrajectory();

	return(0);
}