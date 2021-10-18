//give every supported system type a unique numeric identifier
//this is to be able to selectively compile only the functions that are needed, or at least avoid all runtime checks
#define PES_X3Y1Z1U1 1
#define PES_X4Y2Z1 2
#define PES_X5Y2Z1 3
#define PES_X6Y2Z1 4
#define PES_X4Y1Z1U1 5
#define PES_X3Y2Z1 6
#define PES_X5Y1Z1U1 7

//turned the nasty preprocessor hack into much nicer C++17 constexpr functions
//TODO: change them to consteval when C++20 comes around
constexpr bool IsSystemBraamsSupported(){
	///returns true if the current system type is supported via the Bastian Braaams PES library, false otherwise
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return true;
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return true;
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return true;
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return true;
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return true;
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y1Z1U1) return true;
	return false;
}

constexpr uint32_t gpscore(){
	///core logic of GetPermSymm()
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return 6;
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return 12;
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return 24;
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return 48;
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return 240;
	if constexpr (PES_SYSTEM_TYPE == PES_X6Y2Z1) return 1440;
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y1Z1U1) return 120;
	return 0;
}
constexpr uint32_t GetPermSymm(){
	///this wrapper is needed because constexpr if does not disarm static asserts in the not-taken branch :(
	///this is not a bug, its due to the current C++ standard being a ****
	constexpr uint32_t PS = gpscore();
	static_assert(PS!=0, "System type not recognized");
	return PS;
}
constexpr uint32_t gnacore(){
	///core logic of GetNatoms()
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return 6;
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return 6;
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return 7;
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return 7;
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return 8;
	if constexpr (PES_SYSTEM_TYPE == PES_X6Y2Z1) return 9;
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y1Z1U1) return 8;
	return 0;
}
constexpr uint32_t GetNatoms(){
	///this wrapper is needed because constexpr if does not disarm static asserts in the not-taken branch :(
	///this is not a bug, its due to the current C++ standard being a ****
	constexpr uint32_t NA = gnacore();
	static_assert(NA!=0, "System type not recognized");
	return NA;
}
