import pickle

"""
Holds the unpickled db
"""
with open("BodyDB.pickle", "rb") as file:
    db: dict[int, tuple[tuple[float, float, float, float, float, float]]] = pickle.load(
        file
    )

BODY_DB = db
