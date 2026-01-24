"""
Database connection management
Clean, structured connection handling for PostgreSQL
"""

import os
import logging
from typing import Generator
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.pool import QueuePool

logger = logging.getLogger(__name__)

# Global engine instance
_engine = None
_SessionLocal = None


def get_db_engine():
    """Get or create the database engine"""
    global _engine
    
    if _engine is None:
        database_url = os.environ.get('DATABASE_URL')
        
        if not database_url:
            raise ValueError("DATABASE_URL environment variable not set")
        
        # Use pg8000 driver (pure Python, works on Replit)
        if database_url.startswith('postgresql://'):
            database_url = database_url.replace('postgresql://', 'postgresql+pg8000://')
        
        # Remove sslmode from URL as pg8000 handles SSL differently
        if '?sslmode=' in database_url:
            base_url = database_url.split('?')[0]
            database_url = base_url
        
        _engine = create_engine(
            database_url,
            poolclass=QueuePool,
            pool_size=5,
            max_overflow=10,
            pool_pre_ping=True,
            echo=False
        )
        logger.info("Database engine created successfully")
    
    return _engine


def get_session_factory():
    """Get or create the session factory"""
    global _SessionLocal
    
    if _SessionLocal is None:
        engine = get_db_engine()
        _SessionLocal = sessionmaker(
            autocommit=False,
            autoflush=False,
            bind=engine
        )
    
    return _SessionLocal


def get_db_session() -> Generator[Session, None, None]:
    """Get a database session (generator for dependency injection)"""
    SessionLocal = get_session_factory()
    session = SessionLocal()
    try:
        yield session
    finally:
        session.close()


@contextmanager
def db_session() -> Generator[Session, None, None]:
    """Context manager for database sessions"""
    SessionLocal = get_session_factory()
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def init_database():
    """Initialize database tables"""
    from .models import Base
    
    engine = get_db_engine()
    Base.metadata.create_all(bind=engine)
    logger.info("Database tables created successfully")
    return engine
